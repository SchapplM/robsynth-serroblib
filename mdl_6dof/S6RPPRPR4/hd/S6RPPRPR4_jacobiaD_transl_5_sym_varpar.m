% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:27:34
% EndTime: 2019-02-26 20:27:34
% DurationCPUTime: 0.08s
% Computational Cost: add. (90->25), mult. (188->36), div. (0->0), fcn. (165->8), ass. (0->19)
t56 = qJ(4) + pkin(10);
t54 = sin(t56);
t55 = cos(t56);
t64 = r_i_i_C(1) * t55 - r_i_i_C(2) * t54 + cos(qJ(4)) * pkin(4);
t72 = qJD(4) * t64;
t65 = sin(qJ(4)) * pkin(4) + r_i_i_C(1) * t54 + r_i_i_C(2) * t55;
t63 = t65 * qJD(4);
t71 = -pkin(1) - pkin(2);
t69 = -r_i_i_C(3) - qJ(5) - pkin(7);
t68 = cos(pkin(9));
t66 = pkin(3) + t64;
t57 = sin(pkin(9));
t60 = sin(qJ(1));
t62 = cos(qJ(1));
t50 = t62 * t57 - t60 * t68;
t51 = t60 * t57 + t62 * t68;
t49 = t50 * qJD(1);
t48 = t51 * qJD(1);
t1 = [qJD(2) * t62 - t51 * qJD(5) + t69 * t49 - t66 * t48 - t50 * t63 + (-qJ(2) * t60 + t71 * t62) * qJD(1), qJD(1) * t62, 0, -t65 * t49 - t51 * t72, -t48, 0; t60 * qJD(2) + t50 * qJD(5) + t69 * t48 + t66 * t49 - t51 * t63 + (qJ(2) * t62 + t71 * t60) * qJD(1), qJD(1) * t60, 0, -t65 * t48 + t50 * t72, t49, 0; 0, 0, 0, t63, 0, 0;];
JaD_transl  = t1;
