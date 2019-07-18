% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RRPR2
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
%
% Output:
% JaD_transl [3x4]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4RRPR2_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_jacobiaD_transl_4_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_jacobiaD_transl_4_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RRPR2_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_jacobiaD_transl_4_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:41
% EndTime: 2019-07-18 18:16:41
% DurationCPUTime: 0.04s
% Computational Cost: add. (164->19), mult. (142->27), div. (0->0), fcn. (116->6), ass. (0->19)
t103 = qJ(1) + qJ(2);
t100 = sin(t103);
t101 = cos(t103);
t102 = qJD(1) + qJD(2);
t104 = sin(qJ(4));
t105 = cos(qJ(4));
t108 = qJD(4) * t105;
t109 = qJD(4) * t104;
t91 = -t100 * t109 - t101 * t108 + (t100 * t104 + t101 * t105) * t102;
t110 = t102 * t101;
t111 = t102 * t100;
t92 = t100 * t108 - t101 * t109 + t104 * t110 - t105 * t111;
t115 = t91 * r_i_i_C(1) - t92 * r_i_i_C(2);
t114 = -t92 * r_i_i_C(1) - t91 * r_i_i_C(2);
t113 = -pkin(2) - pkin(3);
t112 = pkin(1) * qJD(1);
t107 = qJ(3) * t110 + t100 * qJD(3) + t113 * t111 - t114;
t106 = t101 * qJD(3) + (-qJ(3) * t100 + t113 * t101) * t102 - t115;
t1 = [-cos(qJ(1)) * t112 + t106, t106, t110, t115; -sin(qJ(1)) * t112 + t107, t107, t111, t114; 0, 0, 0, 0;];
JaD_transl  = t1;
