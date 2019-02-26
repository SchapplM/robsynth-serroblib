% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR15_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR15_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:15
% EndTime: 2019-02-26 22:24:15
% DurationCPUTime: 0.37s
% Computational Cost: add. (589->93), mult. (1627->159), div. (0->0), fcn. (2127->14), ass. (0->63)
t88 = pkin(4) + pkin(10);
t53 = sin(qJ(2));
t54 = sin(qJ(1));
t58 = cos(qJ(2));
t59 = cos(qJ(1));
t71 = cos(pkin(6));
t65 = t59 * t71;
t39 = t53 * t65 + t54 * t58;
t52 = sin(qJ(3));
t57 = cos(qJ(3));
t38 = t54 * t53 - t58 * t65;
t49 = cos(pkin(7));
t47 = sin(pkin(7));
t48 = sin(pkin(6));
t78 = t48 * t59;
t68 = t47 * t78;
t62 = t38 * t49 + t68;
t87 = -t39 * t57 + t62 * t52;
t51 = sin(qJ(5));
t56 = cos(qJ(5));
t50 = sin(qJ(6));
t55 = cos(qJ(6));
t63 = t55 * r_i_i_C(1) - t50 * r_i_i_C(2) + pkin(5);
t85 = r_i_i_C(3) + pkin(12);
t60 = t63 * t51 - t85 * t56 + qJ(4);
t86 = pkin(3) + pkin(11);
t84 = t39 * t52;
t82 = t47 * t48;
t81 = t47 * t51;
t80 = t47 * t56;
t79 = t48 * t54;
t77 = t49 * t52;
t76 = t49 * t57;
t75 = t52 * t53;
t74 = t52 * t58;
t73 = t53 * t57;
t72 = t57 * t58;
t70 = t53 * t82;
t69 = t47 * t79;
t67 = t88 * t47;
t66 = t54 * t71;
t64 = t71 * t47;
t30 = -t38 * t47 + t49 * t78;
t40 = -t59 * t53 - t58 * t66;
t32 = -t40 * t47 + t49 * t79;
t61 = t50 * r_i_i_C(1) + t55 * r_i_i_C(2) + t86;
t41 = -t53 * t66 + t59 * t58;
t37 = t71 * t49 - t58 * t82;
t35 = (t49 * t73 + t74) * t48;
t29 = t52 * t64 + (t49 * t74 + t73) * t48;
t28 = -t57 * t64 + (-t49 * t72 + t75) * t48;
t24 = t40 * t52 + t41 * t76;
t22 = -t38 * t52 + t39 * t76;
t21 = t41 * t57 + (t40 * t49 + t69) * t52;
t20 = -t40 * t76 + t41 * t52 - t57 * t69;
t16 = t38 * t76 + t57 * t68 + t84;
t15 = t28 * t51 + t37 * t56;
t8 = t20 * t51 + t32 * t56;
t7 = -t20 * t56 + t32 * t51;
t6 = t16 * t51 - t30 * t56;
t2 = t21 * t50 + t8 * t55;
t1 = t21 * t55 - t8 * t50;
t3 = [pkin(9) * t78 - t54 * pkin(1) - t39 * pkin(2) + t61 * t87 + t60 * (-t62 * t57 - t84) + (t85 * t51 + t63 * t56 + t88) * t30, t40 * pkin(2) + t24 * qJ(4) + t41 * t67 + t63 * (t24 * t51 + t41 * t80) + t61 * (t40 * t57 - t41 * t77) - t85 * (t24 * t56 - t41 * t81) -t61 * t20 + t60 * t21, t20, -t63 * t7 + t85 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t59 * pkin(1) + t41 * pkin(2) + t8 * pkin(5) + pkin(9) * t79 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t20 * qJ(4) + t86 * t21 + t88 * t32 + t85 * t7, -t38 * pkin(2) + t22 * qJ(4) - t85 * (t22 * t56 - t39 * t81) + t39 * t67 + t63 * (t22 * t51 + t39 * t80) + t61 * (-t38 * t57 - t39 * t77) -t61 * t16 - t60 * t87, t16, t85 * t6 + t63 * (t16 * t56 + t30 * t51) (-t6 * t50 - t55 * t87) * r_i_i_C(1) + (t50 * t87 - t6 * t55) * r_i_i_C(2); 0, t35 * qJ(4) + t63 * (t35 * t51 + t56 * t70) + t85 * (-t35 * t56 + t51 * t70) + (t61 * (-t49 * t75 + t72) + t58 * pkin(2) + t53 * t67) * t48, -t61 * t28 + t60 * t29, t28, t85 * t15 + t63 * (t28 * t56 - t37 * t51) (-t15 * t50 + t29 * t55) * r_i_i_C(1) + (-t15 * t55 - t29 * t50) * r_i_i_C(2);];
Ja_transl  = t3;
