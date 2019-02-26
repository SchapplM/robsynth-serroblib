% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP6_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP6_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:49
% EndTime: 2019-02-26 21:48:50
% DurationCPUTime: 0.23s
% Computational Cost: add. (288->57), mult. (727->94), div. (0->0), fcn. (957->12), ass. (0->42)
t30 = sin(pkin(11));
t35 = sin(qJ(2));
t39 = cos(qJ(2));
t49 = cos(pkin(11));
t24 = -t30 * t39 - t35 * t49;
t36 = sin(qJ(1));
t40 = cos(qJ(1));
t32 = cos(pkin(6));
t43 = -t30 * t35 + t39 * t49;
t42 = t43 * t32;
t10 = t36 * t24 + t40 * t42;
t33 = sin(qJ(5));
t37 = cos(qJ(5));
t21 = t24 * t32;
t11 = -t21 * t40 + t36 * t43;
t34 = sin(qJ(4));
t38 = cos(qJ(4));
t31 = sin(pkin(6));
t50 = t40 * t31;
t4 = t11 * t38 - t34 * t50;
t59 = t10 * t37 + t33 * t4;
t58 = t10 * t33 - t37 * t4;
t46 = r_i_i_C(1) * t37 - r_i_i_C(2) * t33 + pkin(4);
t57 = r_i_i_C(3) + pkin(10);
t41 = t34 * t57 + t38 * t46 + pkin(3);
t56 = t39 * pkin(2);
t53 = t32 * t39;
t51 = t36 * t31;
t47 = -t21 * t36 - t40 * t43;
t45 = r_i_i_C(1) * t33 + r_i_i_C(2) * t37 + pkin(9);
t44 = -t11 * t34 - t38 * t50;
t29 = pkin(1) + t56;
t22 = t32 * t35 * pkin(2) + (-pkin(8) - qJ(3)) * t31;
t20 = t24 * t31;
t19 = t43 * t31;
t16 = -t20 * t38 + t32 * t34;
t13 = t24 * t40 - t36 * t42;
t8 = t34 * t51 - t38 * t47;
t7 = -t34 * t47 - t38 * t51;
t2 = -t13 * t33 + t37 * t8;
t1 = -t13 * t37 - t33 * t8;
t3 = [-t11 * pkin(3) - pkin(4) * t4 + t10 * pkin(9) + r_i_i_C(1) * t58 + r_i_i_C(2) * t59 - t40 * t22 - t36 * t29 + t57 * t44 (-t35 * t40 - t36 * t53) * pkin(2) - t45 * t47 + t41 * t13, t51, -t46 * t7 + t57 * t8, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; -pkin(3) * t47 + t8 * pkin(4) - t13 * pkin(9) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t36 * t22 + t40 * t29 + t57 * t7 (-t35 * t36 + t40 * t53) * pkin(2) + t45 * t11 + t41 * t10, -t50, t4 * t57 + t44 * t46, -r_i_i_C(1) * t59 + r_i_i_C(2) * t58, 0; 0, t19 * t41 - t20 * t45 + t31 * t56, t32, t57 * t16 + t46 * (t20 * t34 + t32 * t38) (-t16 * t33 - t19 * t37) * r_i_i_C(1) + (-t16 * t37 + t19 * t33) * r_i_i_C(2), 0;];
Ja_transl  = t3;
