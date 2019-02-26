% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:29:33
% EndTime: 2019-02-26 21:29:33
% DurationCPUTime: 0.24s
% Computational Cost: add. (384->64), mult. (758->96), div. (0->0), fcn. (998->14), ass. (0->45)
t35 = sin(pkin(11));
t40 = sin(qJ(2));
t43 = cos(qJ(2));
t54 = cos(pkin(11));
t24 = -t43 * t35 - t40 * t54;
t41 = sin(qJ(1));
t44 = cos(qJ(1));
t37 = cos(pkin(6));
t47 = -t40 * t35 + t43 * t54;
t46 = t47 * t37;
t10 = t41 * t24 + t44 * t46;
t39 = sin(qJ(6));
t21 = t24 * t37;
t11 = -t44 * t21 + t41 * t47;
t33 = pkin(12) + qJ(5);
t31 = sin(t33);
t32 = cos(t33);
t36 = sin(pkin(6));
t55 = t44 * t36;
t4 = t11 * t32 - t31 * t55;
t42 = cos(qJ(6));
t64 = t10 * t42 + t4 * t39;
t63 = t10 * t39 - t4 * t42;
t29 = cos(pkin(12)) * pkin(4) + pkin(3);
t50 = t42 * r_i_i_C(1) - t39 * r_i_i_C(2) + pkin(5);
t62 = pkin(10) + r_i_i_C(3);
t45 = t62 * t31 + t50 * t32 + t29;
t61 = t43 * pkin(2);
t58 = t37 * t43;
t56 = t41 * t36;
t52 = -t37 * t40 * pkin(2) + (sin(pkin(12)) * pkin(4) + pkin(8) + qJ(3)) * t36;
t51 = -t41 * t21 - t44 * t47;
t38 = -pkin(9) - qJ(4);
t49 = t39 * r_i_i_C(1) + t42 * r_i_i_C(2) - t38;
t48 = -t11 * t31 - t32 * t55;
t30 = pkin(1) + t61;
t20 = t24 * t36;
t19 = t47 * t36;
t16 = -t20 * t32 + t37 * t31;
t13 = t44 * t24 - t41 * t46;
t8 = t31 * t56 - t32 * t51;
t7 = -t31 * t51 - t32 * t56;
t2 = -t13 * t39 + t8 * t42;
t1 = -t13 * t42 - t8 * t39;
t3 = [-t4 * pkin(5) + t63 * r_i_i_C(1) + t64 * r_i_i_C(2) - t10 * t38 - t11 * t29 - t41 * t30 + t52 * t44 + t62 * t48 (-t44 * t40 - t41 * t58) * pkin(2) - t49 * t51 + t45 * t13, t56, -t13, -t50 * t7 + t62 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t8 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t13 * t38 - t29 * t51 + t44 * t30 + t52 * t41 + t62 * t7 (-t41 * t40 + t44 * t58) * pkin(2) + t49 * t11 + t45 * t10, -t55, -t10, t62 * t4 + t50 * t48, -t64 * r_i_i_C(1) + t63 * r_i_i_C(2); 0, t45 * t19 - t49 * t20 + t36 * t61, t37, -t19, t62 * t16 + t50 * (t20 * t31 + t37 * t32) (-t16 * t39 - t19 * t42) * r_i_i_C(1) + (-t16 * t42 + t19 * t39) * r_i_i_C(2);];
Ja_transl  = t3;
