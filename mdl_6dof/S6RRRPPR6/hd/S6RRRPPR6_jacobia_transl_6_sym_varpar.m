% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:06:22
% EndTime: 2019-02-26 22:06:23
% DurationCPUTime: 0.19s
% Computational Cost: add. (302->55), mult. (498->89), div. (0->0), fcn. (629->12), ass. (0->38)
t29 = cos(qJ(3));
t18 = t29 * pkin(3) + pkin(2);
t21 = qJ(3) + pkin(11);
t19 = sin(t21);
t20 = cos(t21);
t24 = sin(qJ(6));
t28 = cos(qJ(6));
t35 = t24 * r_i_i_C(1) + t28 * r_i_i_C(2) + qJ(5);
t39 = pkin(4) + pkin(10) + r_i_i_C(3);
t46 = t35 * t19 + t39 * t20 + t18;
t45 = pkin(5) + qJ(4) + pkin(9);
t22 = sin(pkin(6));
t26 = sin(qJ(2));
t44 = t22 * t26;
t27 = sin(qJ(1));
t43 = t22 * t27;
t30 = cos(qJ(2));
t42 = t22 * t30;
t31 = cos(qJ(1));
t41 = t22 * t31;
t40 = cos(pkin(6));
t38 = t27 * t40;
t37 = t31 * t40;
t25 = sin(qJ(3));
t36 = t22 * (pkin(3) * t25 + pkin(8));
t12 = t26 * t37 + t27 * t30;
t3 = t12 * t19 + t20 * t41;
t34 = -t12 * t20 + t19 * t41;
t33 = t28 * r_i_i_C(1) - t24 * r_i_i_C(2) + t45;
t14 = -t26 * t38 + t31 * t30;
t13 = t31 * t26 + t30 * t38;
t11 = t27 * t26 - t30 * t37;
t9 = t19 * t44 - t40 * t20;
t8 = t14 * t20 + t19 * t43;
t7 = t14 * t19 - t20 * t43;
t2 = t13 * t28 + t7 * t24;
t1 = -t13 * t24 + t7 * t28;
t4 = [-t27 * pkin(1) - t33 * t11 - t12 * t18 - t35 * t3 + t31 * t36 + t39 * t34, -t13 * t46 + t33 * t14 (-t14 * t25 + t29 * t43) * pkin(3) + t35 * t8 - t39 * t7, t13, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t31 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * qJ(5) + t45 * t13 + t14 * t18 + t27 * t36 + t39 * t8, -t11 * t46 + t33 * t12 (-t12 * t25 - t29 * t41) * pkin(3) - t35 * t34 - t39 * t3, t11, t3 (-t11 * t24 + t3 * t28) * r_i_i_C(1) + (-t11 * t28 - t3 * t24) * r_i_i_C(2); 0 (t33 * t26 + t46 * t30) * t22 (-t25 * t44 + t40 * t29) * pkin(3) - t39 * t9 + t35 * (t40 * t19 + t20 * t44) -t42, t9 (t24 * t42 + t9 * t28) * r_i_i_C(1) + (-t9 * t24 + t28 * t42) * r_i_i_C(2);];
Ja_transl  = t4;
