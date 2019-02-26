% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:16:26
% EndTime: 2019-02-26 21:16:26
% DurationCPUTime: 0.13s
% Computational Cost: add. (269->40), mult. (176->50), div. (0->0), fcn. (115->9), ass. (0->38)
t56 = pkin(11) + qJ(3);
t51 = sin(t56);
t57 = qJD(3) + qJD(4);
t53 = qJD(5) + t57;
t54 = qJ(4) + t56;
t50 = qJ(5) + t54;
t47 = cos(t50);
t79 = r_i_i_C(2) * t47;
t46 = sin(t50);
t81 = r_i_i_C(1) * t46;
t64 = t79 + t81;
t62 = t64 * t53;
t48 = sin(t54);
t82 = pkin(4) * t48;
t60 = -t57 * t82 - t62;
t75 = pkin(3) * qJD(3);
t84 = -t51 * t75 + t60;
t77 = t47 * t53;
t70 = r_i_i_C(1) * t77;
t49 = cos(t54);
t76 = t49 * t57;
t83 = -pkin(4) * t76 - t70;
t80 = r_i_i_C(2) * t46;
t78 = r_i_i_C(3) + pkin(9) + pkin(8) + pkin(7) + qJ(2);
t58 = sin(qJ(1));
t74 = qJD(1) * t58;
t59 = cos(qJ(1));
t73 = qJD(1) * t59;
t69 = t53 * t80;
t66 = qJD(1) * t79;
t68 = t58 * t66 + t59 * t69 + t74 * t81;
t52 = cos(t56);
t65 = -t52 * t75 + t83;
t63 = -r_i_i_C(1) * t47 - pkin(4) * t49 - pkin(3) * t52 - cos(pkin(11)) * pkin(2) - pkin(1) + t80;
t61 = -t59 * t70 + t68;
t43 = -pkin(3) * t51 - t82;
t41 = t58 * t69;
t1 = [t59 * qJD(2) - t84 * t58 + (-t78 * t58 + t63 * t59) * qJD(1), t73, -t43 * t74 + t65 * t59 + t68 (t48 * t74 - t59 * t76) * pkin(4) + t61, t61, 0; t58 * qJD(2) + t84 * t59 + (t63 * t58 + t78 * t59) * qJD(1), t74, t41 + t65 * t58 + (t43 - t64) * t73, t41 + t83 * t58 + (-t64 - t82) * t73, -t59 * t66 + t41 + (-t46 * t73 - t58 * t77) * r_i_i_C(1), 0; 0, 0, t84, t60, -t62, 0;];
JaD_transl  = t1;
