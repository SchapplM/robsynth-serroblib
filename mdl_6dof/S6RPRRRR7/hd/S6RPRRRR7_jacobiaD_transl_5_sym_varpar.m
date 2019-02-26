% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR7_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR7_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:18:11
% EndTime: 2019-02-26 21:18:11
% DurationCPUTime: 0.11s
% Computational Cost: add. (206->39), mult. (182->51), div. (0->0), fcn. (117->8), ass. (0->40)
t56 = qJ(3) + qJ(4);
t53 = qJ(5) + t56;
t48 = sin(t53);
t82 = r_i_i_C(2) * t48;
t52 = cos(t56);
t84 = pkin(4) * t52;
t85 = -t82 + t84;
t49 = cos(t53);
t83 = r_i_i_C(1) * t49;
t54 = qJD(3) + qJD(4);
t50 = qJD(5) + t54;
t81 = t48 * t50;
t80 = t49 * t50;
t51 = sin(t56);
t79 = t51 * t54;
t78 = pkin(3) * qJD(3);
t58 = sin(qJ(1));
t77 = qJD(1) * t58;
t60 = cos(qJ(1));
t76 = qJD(1) * t60;
t75 = -pkin(1) - r_i_i_C(3) - pkin(9) - pkin(8) - pkin(7);
t74 = pkin(4) * t79;
t73 = t54 * t84;
t72 = r_i_i_C(1) * t81;
t70 = qJD(1) * t83;
t71 = t58 * t70 + (t80 * r_i_i_C(2) + t72) * t60;
t59 = cos(qJ(3));
t69 = t59 * t78;
t68 = -r_i_i_C(1) * t80 + r_i_i_C(2) * t81;
t67 = -r_i_i_C(1) * t48 - r_i_i_C(2) * t49;
t66 = qJD(1) * (t59 * pkin(3) + t85);
t65 = t67 * t50;
t57 = sin(qJ(3));
t64 = t57 * pkin(3) + pkin(4) * t51 + qJ(2) - t67;
t63 = -t77 * t82 + t71;
t62 = t68 - t73;
t61 = qJD(2) + t69 + t73 + (-t82 + t83) * t50;
t45 = t60 * t70;
t39 = -t57 * t78 - t74;
t1 = [t61 * t60 + (-t64 * t58 + t75 * t60) * qJD(1), t76, t45 + t60 * t66 + (t39 + t65) * t58, t45 + t85 * t76 + (t65 - t74) * t58, -t58 * t72 + t45 + (-t48 * t76 - t58 * t80) * r_i_i_C(2), 0; t61 * t58 + (t75 * t58 + t64 * t60) * qJD(1), t77, -t60 * t39 + t58 * t66 + t71 (t52 * t77 + t60 * t79) * pkin(4) + t63, t63, 0; 0, 0, t62 - t69, t62, t68, 0;];
JaD_transl  = t1;
