% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP9
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRP9_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP9_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:44:25
% EndTime: 2019-02-26 22:44:25
% DurationCPUTime: 0.20s
% Computational Cost: add. (331->58), mult. (576->91), div. (0->0), fcn. (723->12), ass. (0->43)
t35 = sin(qJ(3));
t38 = cos(qJ(3));
t33 = qJ(4) + qJ(5);
t30 = cos(t33);
t23 = pkin(5) * t30 + cos(qJ(4)) * pkin(4);
t21 = pkin(3) + t23;
t29 = sin(t33);
t43 = r_i_i_C(1) * t30 - r_i_i_C(2) * t29 + t21;
t53 = r_i_i_C(3) + qJ(6) + pkin(11) + pkin(10);
t58 = t53 * t35 + t43 * t38 + pkin(2);
t36 = sin(qJ(2));
t37 = sin(qJ(1));
t39 = cos(qJ(2));
t48 = cos(pkin(6));
t52 = cos(qJ(1));
t45 = t48 * t52;
t18 = t36 * t45 + t37 * t39;
t34 = sin(pkin(6));
t47 = t34 * t52;
t10 = t18 * t38 - t35 * t47;
t17 = t37 * t36 - t39 * t45;
t44 = -t10 * t29 + t17 * t30;
t57 = t44 * r_i_i_C(1) + (-t10 * t30 - t17 * t29) * r_i_i_C(2);
t46 = t37 * t48;
t20 = -t36 * t46 + t52 * t39;
t51 = t34 * t37;
t14 = t20 * t38 + t35 * t51;
t19 = t52 * t36 + t39 * t46;
t5 = -t14 * t29 + t19 * t30;
t6 = t14 * t30 + t19 * t29;
t56 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t50 = t34 * t38;
t16 = t48 * t35 + t36 * t50;
t49 = t34 * t39;
t42 = -t16 * t29 - t30 * t49;
t55 = t42 * r_i_i_C(1) + (-t16 * t30 + t29 * t49) * r_i_i_C(2);
t22 = pkin(5) * t29 + sin(qJ(4)) * pkin(4);
t54 = pkin(9) + t22;
t41 = t29 * r_i_i_C(1) + t30 * r_i_i_C(2) + t54;
t9 = t18 * t35 + t38 * t47;
t15 = t34 * t36 * t35 - t48 * t38;
t13 = t20 * t35 - t37 * t50;
t1 = [-t37 * pkin(1) - t18 * pkin(2) + pkin(8) * t47 - t43 * t10 - t41 * t17 - t53 * t9, -t19 * t58 + t41 * t20, -t43 * t13 + t53 * t14, -t14 * t22 + t19 * t23 + t56, t5 * pkin(5) + t56, t13; t52 * pkin(1) + t20 * pkin(2) + pkin(8) * t51 + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t53 * t13 + t14 * t21 + t54 * t19, -t17 * t58 + t41 * t18, t53 * t10 - t43 * t9, -t10 * t22 + t17 * t23 + t57, t44 * pkin(5) + t57, t9; 0 (t41 * t36 + t58 * t39) * t34, -t43 * t15 + t53 * t16, -t16 * t22 - t23 * t49 + t55, t42 * pkin(5) + t55, t15;];
Ja_transl  = t1;
