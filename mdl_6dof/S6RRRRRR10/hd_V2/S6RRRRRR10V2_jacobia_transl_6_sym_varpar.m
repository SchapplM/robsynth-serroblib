% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR10V2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10V2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobia_transl_6_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:56:31
% EndTime: 2019-04-11 14:56:32
% DurationCPUTime: 0.30s
% Computational Cost: add. (342->62), mult. (503->112), div. (0->0), fcn. (621->12), ass. (0->56)
t42 = qJ(2) + qJ(3);
t40 = cos(t42);
t39 = sin(t42);
t49 = cos(qJ(5));
t44 = sin(qJ(5));
t50 = cos(qJ(4));
t73 = t44 * t50;
t54 = t39 * t73 + t40 * t49;
t81 = r_i_i_C(3) + pkin(6);
t86 = pkin(5) * t40 - t54 * t81;
t43 = sin(qJ(6));
t48 = cos(qJ(6));
t60 = r_i_i_C(1) * t48 - r_i_i_C(2) * t43;
t45 = sin(qJ(4));
t51 = cos(qJ(1));
t67 = t51 * t45;
t47 = sin(qJ(1));
t69 = t47 * t50;
t27 = t40 * t69 - t67;
t77 = t39 * t47;
t10 = t27 * t49 + t44 * t77;
t66 = t51 * t50;
t70 = t47 * t45;
t26 = t40 * t70 + t66;
t84 = t10 * t43 - t26 * t48;
t83 = -t10 * t48 - t26 * t43;
t82 = pkin(3) * t40 + pkin(5) * t39;
t76 = t39 * t49;
t75 = t39 * t51;
t74 = t43 * t45;
t72 = t45 * t48;
t71 = t45 * t49;
t68 = t49 * t50;
t65 = t39 * t74;
t64 = t39 * t72;
t63 = t39 * t67;
t62 = t81 * t44;
t61 = -sin(qJ(2)) * pkin(2) - pkin(3) * t39;
t59 = r_i_i_C(1) * t43 + r_i_i_C(2) * t48;
t41 = cos(qJ(2)) * pkin(2);
t58 = t41 + pkin(1) + t82;
t23 = t39 * t68 - t40 * t44;
t57 = (-r_i_i_C(1) * t65 - r_i_i_C(2) * t64 - t23 * t60 + t86) * t47;
t20 = t23 * t51;
t56 = (t20 * t43 - t48 * t63) * r_i_i_C(2) + (-t20 * t48 - t43 * t63) * r_i_i_C(1) + t86 * t51;
t55 = -t27 * t44 + t47 * t76;
t25 = t39 * t44 + t40 * t68;
t53 = (-t25 * t43 + t40 * t72) * r_i_i_C(2) + (t25 * t48 + t40 * t74) * r_i_i_C(1) + t82 + t81 * (t40 * t73 - t76);
t52 = -t49 * t60 - t62;
t29 = t40 * t66 + t70;
t28 = t40 * t67 - t69;
t14 = t29 * t49 + t44 * t75;
t13 = t29 * t44 - t49 * t75;
t2 = t14 * t48 + t28 * t43;
t1 = -t14 * t43 + t28 * t48;
t3 = [r_i_i_C(1) * t83 + r_i_i_C(2) * t84 - t58 * t47 + t81 * t55, t51 * t61 + t56, -pkin(3) * t75 + t56, t28 * t52 + t29 * t59, -t13 * t60 + t14 * t81, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t13 * t81 + t51 * t58, t47 * t61 + t57, -pkin(3) * t77 + t57, t26 * t52 + t27 * t59, t10 * t81 + t55 * t60, -r_i_i_C(1) * t84 + r_i_i_C(2) * t83; 0, t41 + t53, t53 ((t43 * t50 - t48 * t71) * r_i_i_C(1) + (t43 * t71 + t48 * t50) * r_i_i_C(2) - t45 * t62) * t39, t23 * t81 - t54 * t60 (-t23 * t43 + t64) * r_i_i_C(1) + (-t23 * t48 - t65) * r_i_i_C(2);];
Ja_transl  = t3;
