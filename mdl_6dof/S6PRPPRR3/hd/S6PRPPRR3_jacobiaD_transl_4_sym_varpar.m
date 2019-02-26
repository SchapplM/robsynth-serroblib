% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPPRR3_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR3_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobiaD_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:50
% EndTime: 2019-02-26 19:45:50
% DurationCPUTime: 0.04s
% Computational Cost: add. (39->16), mult. (133->30), div. (0->0), fcn. (120->8), ass. (0->20)
t103 = cos(qJ(2));
t97 = sin(pkin(10));
t112 = t97 * t103;
t100 = cos(pkin(10));
t111 = t100 * t103;
t101 = cos(pkin(6));
t102 = sin(qJ(2));
t110 = t101 * t102;
t109 = qJD(2) * t102;
t108 = t97 * t109;
t107 = qJD(2) * t111;
t96 = sin(pkin(11));
t99 = cos(pkin(11));
t106 = t96 * r_i_i_C(1) + t99 * r_i_i_C(2) + qJ(3);
t105 = -t99 * r_i_i_C(1) + t96 * r_i_i_C(2) - pkin(2) - pkin(3);
t104 = t100 * t110 + t112;
t98 = sin(pkin(6));
t93 = -t101 * t108 + t107;
t91 = t104 * qJD(2);
t1 = [0 -(t97 * t110 - t111) * qJD(3) - t106 * (t100 * t102 + t101 * t112) * qJD(2) + t105 * t93, t93, 0, 0, 0; 0, t104 * qJD(3) - t106 * (-t101 * t107 + t108) + t105 * t91, t91, 0, 0, 0; 0 (t102 * qJD(3) + (t105 * t102 + t106 * t103) * qJD(2)) * t98, t98 * t109, 0, 0, 0;];
JaD_transl  = t1;
