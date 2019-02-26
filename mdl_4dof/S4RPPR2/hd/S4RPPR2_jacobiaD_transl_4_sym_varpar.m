% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RPPR2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
%
% Output:
% JaD_transl [3x4]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4RPPR2_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_jacobiaD_transl_4_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_jacobiaD_transl_4_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPPR2_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_jacobiaD_transl_4_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:31:23
% EndTime: 2019-02-26 19:31:23
% DurationCPUTime: 0.04s
% Computational Cost: add. (72->18), mult. (96->25), div. (0->0), fcn. (80->6), ass. (0->16)
t98 = -pkin(1) - cos(pkin(6)) * pkin(3) - pkin(2);
t89 = sin(qJ(1));
t97 = qJD(1) * t89;
t90 = cos(qJ(1));
t96 = qJD(1) * t90;
t95 = qJD(4) * t89;
t94 = qJD(4) * t90;
t93 = pkin(3) * sin(pkin(6)) + qJ(2);
t87 = pkin(6) + qJ(4);
t85 = sin(t87);
t86 = cos(t87);
t78 = -t85 * t95 - t86 * t94 + (t85 * t89 + t86 * t90) * qJD(1);
t79 = (t95 - t97) * t86 + (-t94 + t96) * t85;
t92 = -t79 * r_i_i_C(1) - t78 * r_i_i_C(2);
t91 = t78 * r_i_i_C(1) - t79 * r_i_i_C(2);
t1 = [t90 * qJD(2) + (-t93 * t89 + t98 * t90) * qJD(1) - t91, t96, 0, t91; t89 * qJD(2) + (t98 * t89 + t93 * t90) * qJD(1) - t92, t97, 0, t92; 0, 0, 0, 0;];
JaD_transl  = t1;
