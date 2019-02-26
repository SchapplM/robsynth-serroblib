% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:05:50
% EndTime: 2019-02-26 22:05:50
% DurationCPUTime: 0.20s
% Computational Cost: add. (324->58), mult. (466->91), div. (0->0), fcn. (586->14), ass. (0->41)
t34 = cos(qJ(3));
t20 = t34 * pkin(3) + pkin(2);
t26 = qJ(3) + pkin(11);
t22 = sin(t26);
t24 = cos(t26);
t19 = cos(pkin(12)) * pkin(5) + pkin(4);
t25 = pkin(12) + qJ(6);
t21 = sin(t25);
t23 = cos(t25);
t39 = t23 * r_i_i_C(1) - t21 * r_i_i_C(2) + t19;
t49 = r_i_i_C(3) + pkin(10) + qJ(5);
t50 = t49 * t22 + t39 * t24 + t20;
t28 = sin(pkin(6));
t32 = sin(qJ(2));
t48 = t28 * t32;
t33 = sin(qJ(1));
t47 = t28 * t33;
t35 = cos(qJ(2));
t46 = t28 * t35;
t36 = cos(qJ(1));
t45 = t28 * t36;
t44 = cos(pkin(6));
t43 = sin(pkin(12)) * pkin(5) + qJ(4) + pkin(9);
t41 = t36 * t44;
t12 = t32 * t41 + t33 * t35;
t4 = t12 * t24 - t22 * t45;
t42 = t33 * t44;
t31 = sin(qJ(3));
t40 = t28 * (pkin(3) * t31 + pkin(8));
t3 = t12 * t22 + t24 * t45;
t38 = t21 * r_i_i_C(1) + t23 * r_i_i_C(2) + t43;
t14 = -t32 * t42 + t36 * t35;
t13 = t36 * t32 + t35 * t42;
t11 = t33 * t32 - t35 * t41;
t10 = t44 * t22 + t24 * t48;
t9 = t22 * t48 - t44 * t24;
t8 = t14 * t24 + t22 * t47;
t7 = t14 * t22 - t24 * t47;
t2 = t13 * t21 + t8 * t23;
t1 = t13 * t23 - t8 * t21;
t5 = [-t33 * pkin(1) - t38 * t11 - t12 * t20 - t49 * t3 + t36 * t40 - t39 * t4, -t13 * t50 + t38 * t14, t49 * t8 + (-t14 * t31 + t34 * t47) * pkin(3) - t39 * t7, t13, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t36 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t43 * t13 + t14 * t20 + t8 * t19 + t33 * t40 + t49 * t7, -t11 * t50 + t38 * t12, t49 * t4 + (-t12 * t31 - t34 * t45) * pkin(3) - t39 * t3, t11, t3 (t11 * t23 - t4 * t21) * r_i_i_C(1) + (-t11 * t21 - t4 * t23) * r_i_i_C(2); 0 (t38 * t32 + t50 * t35) * t28, t49 * t10 + (-t31 * t48 + t44 * t34) * pkin(3) - t39 * t9, -t46, t9 (-t10 * t21 - t23 * t46) * r_i_i_C(1) + (-t10 * t23 + t21 * t46) * r_i_i_C(2);];
Ja_transl  = t5;
