% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRP6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:18:01
% EndTime: 2019-02-26 20:18:02
% DurationCPUTime: 0.34s
% Computational Cost: add. (646->96), mult. (1814->175), div. (0->0), fcn. (2382->14), ass. (0->65)
t57 = sin(pkin(7));
t63 = sin(qJ(2));
t66 = cos(qJ(2));
t59 = cos(pkin(12));
t84 = cos(pkin(6));
t78 = t59 * t84;
t82 = sin(pkin(12));
t69 = t82 * t63 - t66 * t78;
t83 = cos(pkin(7));
t58 = sin(pkin(6));
t90 = t58 * t59;
t98 = t57 * t90 + t69 * t83;
t73 = t84 * t82;
t70 = t59 * t63 + t66 * t73;
t79 = t58 * t82;
t97 = -t57 * t79 + t70 * t83;
t96 = r_i_i_C(1) + pkin(5);
t95 = r_i_i_C(2) + pkin(11);
t94 = pkin(9) * t57;
t93 = cos(qJ(3));
t61 = sin(qJ(4));
t92 = t57 * t61;
t65 = cos(qJ(4));
t91 = t57 * t65;
t89 = t58 * t63;
t88 = t58 * t66;
t60 = sin(qJ(5));
t87 = t60 * t65;
t64 = cos(qJ(5));
t86 = t64 * t65;
t85 = r_i_i_C(3) + qJ(6);
t81 = t57 * t89;
t62 = sin(qJ(3));
t77 = t62 * t83;
t76 = t84 * t57;
t74 = t83 * t93;
t72 = -pkin(4) * t65 - t95 * t61 - pkin(3);
t71 = t85 * t60 + t96 * t64 + pkin(4);
t52 = t59 * t66 - t63 * t73;
t51 = t63 * t78 + t82 * t66;
t50 = -t57 * t88 + t84 * t83;
t49 = (-t63 * t77 + t93 * t66) * t58;
t48 = (t62 * t66 + t63 * t74) * t58;
t45 = t70 * t57 + t83 * t79;
t44 = t69 * t57 - t83 * t90;
t43 = t62 * t76 + (t93 * t63 + t66 * t77) * t58;
t42 = t62 * t89 - t74 * t88 - t93 * t76;
t40 = t49 * t65 + t61 * t81;
t38 = -t52 * t77 - t70 * t93;
t37 = t52 * t74 - t70 * t62;
t36 = -t51 * t77 - t69 * t93;
t35 = t51 * t74 - t69 * t62;
t34 = t43 * t65 + t50 * t61;
t32 = t52 * t93 - t97 * t62;
t31 = t52 * t62 + t97 * t93;
t30 = t51 * t93 - t98 * t62;
t29 = t51 * t62 + t98 * t93;
t24 = t38 * t65 + t52 * t92;
t22 = t36 * t65 + t51 * t92;
t18 = t32 * t65 + t45 * t61;
t16 = t30 * t65 + t44 * t61;
t13 = t34 * t60 - t42 * t64;
t3 = t18 * t60 - t31 * t64;
t1 = t16 * t60 - t29 * t64;
t2 = [0, t24 * pkin(4) + t38 * pkin(3) + t37 * pkin(10) - t70 * pkin(2) + t52 * t94 + t95 * (t38 * t61 - t52 * t91) + t96 * (t24 * t64 + t37 * t60) + t85 * (t24 * t60 - t37 * t64) t32 * pkin(10) + t96 * (-t31 * t86 + t32 * t60) + t85 * (-t31 * t87 - t32 * t64) + t72 * t31, t95 * t18 + t71 * (-t32 * t61 + t45 * t65) t85 * (t18 * t64 + t31 * t60) - t96 * t3, t3; 0, t22 * pkin(4) + t36 * pkin(3) + t35 * pkin(10) - t69 * pkin(2) + t51 * t94 + t85 * (t22 * t60 - t35 * t64) + t95 * (t36 * t61 - t51 * t91) + t96 * (t22 * t64 + t35 * t60) t30 * pkin(10) + t96 * (-t29 * t86 + t30 * t60) + t85 * (-t29 * t87 - t30 * t64) + t72 * t29, t95 * t16 + t71 * (-t30 * t61 + t44 * t65) t85 * (t16 * t64 + t29 * t60) - t96 * t1, t1; 1, t49 * pkin(3) + t40 * pkin(4) + t48 * pkin(10) + (pkin(2) * t66 + t63 * t94) * t58 + t95 * (t49 * t61 - t65 * t81) + t96 * (t40 * t64 + t48 * t60) + t85 * (t40 * t60 - t48 * t64) t43 * pkin(10) + t96 * (-t42 * t86 + t43 * t60) + t85 * (-t42 * t87 - t43 * t64) + t72 * t42, t95 * t34 + t71 * (-t43 * t61 + t50 * t65) t85 * (t34 * t64 + t42 * t60) - t96 * t13, t13;];
Ja_transl  = t2;
