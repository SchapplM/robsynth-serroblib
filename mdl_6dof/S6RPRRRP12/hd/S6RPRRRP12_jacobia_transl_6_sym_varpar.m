% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRP12_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP12_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:14:17
% EndTime: 2019-02-26 21:14:18
% DurationCPUTime: 0.35s
% Computational Cost: add. (659->76), mult. (1819->128), div. (0->0), fcn. (2406->14), ass. (0->55)
t54 = cos(qJ(1));
t70 = sin(pkin(12));
t74 = cos(pkin(6));
t62 = t74 * t70;
t72 = cos(pkin(12));
t79 = sin(qJ(1));
t41 = t54 * t62 + t79 * t72;
t51 = sin(qJ(3));
t80 = cos(qJ(3));
t64 = t74 * t72;
t40 = -t54 * t64 + t79 * t70;
t48 = sin(pkin(6));
t71 = sin(pkin(7));
t66 = t48 * t71;
t73 = cos(pkin(7));
t85 = t40 * t73 + t54 * t66;
t26 = t41 * t80 - t85 * t51;
t67 = t48 * t73;
t36 = t40 * t71 - t54 * t67;
t50 = sin(qJ(4));
t53 = cos(qJ(4));
t14 = t26 * t53 + t36 * t50;
t25 = t41 * t51 + t85 * t80;
t49 = sin(qJ(5));
t52 = cos(qJ(5));
t1 = t14 * t49 - t25 * t52;
t88 = t14 * t52 + t25 * t49;
t84 = -t26 * t50 + t36 * t53;
t57 = t54 * t70 + t79 * t64;
t83 = t57 * t73 - t79 * t66;
t75 = r_i_i_C(3) + qJ(6);
t82 = pkin(5) + r_i_i_C(1);
t58 = t75 * t49 + t82 * t52 + pkin(4);
t81 = pkin(11) + r_i_i_C(2);
t78 = t49 * t53;
t77 = t52 * t53;
t76 = t54 * t48;
t69 = t79 * t48;
t63 = t74 * t71;
t61 = t73 * t72;
t59 = -pkin(4) * t53 - t81 * t50 - pkin(3);
t55 = t57 * t71 + t79 * t67;
t42 = t54 * t72 - t79 * t62;
t39 = -t72 * t66 + t74 * t73;
t34 = t51 * t63 + (t51 * t61 + t80 * t70) * t48;
t33 = -t80 * t63 + (t51 * t70 - t61 * t80) * t48;
t30 = t42 * t80 - t83 * t51;
t29 = t42 * t51 + t83 * t80;
t24 = t34 * t53 + t39 * t50;
t18 = t30 * t53 + t50 * t55;
t17 = t30 * t50 - t53 * t55;
t11 = t24 * t49 - t33 * t52;
t6 = t18 * t52 + t29 * t49;
t5 = t18 * t49 - t29 * t52;
t2 = [-t79 * pkin(1) - t41 * pkin(2) - t26 * pkin(3) - t14 * pkin(4) - t36 * pkin(9) - t25 * pkin(10) + qJ(2) * t76 - t75 * t1 + t81 * t84 - t82 * t88, t69, t30 * pkin(10) + t75 * (-t29 * t78 - t30 * t52) + t82 * (-t29 * t77 + t30 * t49) + t59 * t29, -t58 * t17 + t81 * t18, -t82 * t5 + t75 * t6, t5; t54 * pkin(1) + t42 * pkin(2) + t30 * pkin(3) + t18 * pkin(4) + t55 * pkin(9) + t29 * pkin(10) + qJ(2) * t69 + t81 * t17 + t75 * t5 + t82 * t6, -t76, t26 * pkin(10) + t82 * (-t25 * t77 + t26 * t49) + t75 * (-t25 * t78 - t26 * t52) + t59 * t25, t81 * t14 + t58 * t84, -t82 * t1 + t75 * t88, t1; 0, t74, t34 * pkin(10) + t82 * (-t33 * t77 + t34 * t49) + t75 * (-t33 * t78 - t34 * t52) + t59 * t33, t81 * t24 + t58 * (-t34 * t50 + t39 * t53) t75 * (t24 * t52 + t33 * t49) - t82 * t11, t11;];
Ja_transl  = t2;
