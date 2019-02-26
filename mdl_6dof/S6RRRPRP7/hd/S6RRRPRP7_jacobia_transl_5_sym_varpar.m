% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP7_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP7_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:12:38
% EndTime: 2019-02-26 22:12:38
% DurationCPUTime: 0.19s
% Computational Cost: add. (250->57), mult. (420->91), div. (0->0), fcn. (528->12), ass. (0->40)
t26 = sin(qJ(2));
t27 = sin(qJ(1));
t30 = cos(qJ(2));
t31 = cos(qJ(1));
t39 = cos(pkin(6));
t37 = t31 * t39;
t11 = t27 * t26 - t30 * t37;
t24 = sin(qJ(5));
t28 = cos(qJ(5));
t12 = t26 * t37 + t27 * t30;
t21 = qJ(3) + pkin(11);
t19 = sin(t21);
t20 = cos(t21);
t22 = sin(pkin(6));
t40 = t22 * t31;
t4 = t12 * t20 - t19 * t40;
t49 = -t11 * t28 + t4 * t24;
t48 = -t11 * t24 - t4 * t28;
t29 = cos(qJ(3));
t18 = t29 * pkin(3) + pkin(2);
t35 = t28 * r_i_i_C(1) - t24 * r_i_i_C(2) + pkin(4);
t46 = pkin(10) + r_i_i_C(3);
t47 = t46 * t19 + t35 * t20 + t18;
t43 = t22 * t26;
t42 = t22 * t27;
t41 = t22 * t30;
t38 = t27 * t39;
t25 = sin(qJ(3));
t36 = t22 * (pkin(3) * t25 + pkin(8));
t23 = -qJ(4) - pkin(9);
t34 = t24 * r_i_i_C(1) + t28 * r_i_i_C(2) - t23;
t33 = -t12 * t19 - t20 * t40;
t14 = -t26 * t38 + t31 * t30;
t13 = t31 * t26 + t30 * t38;
t10 = t39 * t19 + t20 * t43;
t8 = t14 * t20 + t19 * t42;
t7 = t14 * t19 - t20 * t42;
t2 = t13 * t24 + t8 * t28;
t1 = t13 * t28 - t8 * t24;
t3 = [-t27 * pkin(1) - t4 * pkin(4) + t48 * r_i_i_C(1) + t49 * r_i_i_C(2) + t11 * t23 - t12 * t18 + t31 * t36 + t46 * t33, -t13 * t47 + t34 * t14, t46 * t8 + (-t14 * t25 + t29 * t42) * pkin(3) - t35 * t7, t13, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t31 * pkin(1) + t8 * pkin(4) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t13 * t23 + t14 * t18 + t27 * t36 + t46 * t7, -t11 * t47 + t34 * t12, t46 * t4 + (-t12 * t25 - t29 * t40) * pkin(3) + t35 * t33, t11, -t49 * r_i_i_C(1) + t48 * r_i_i_C(2), 0; 0 (t34 * t26 + t47 * t30) * t22, t46 * t10 + (-t25 * t43 + t39 * t29) * pkin(3) + t35 * (-t19 * t43 + t39 * t20) -t41 (-t10 * t24 - t28 * t41) * r_i_i_C(1) + (-t10 * t28 + t24 * t41) * r_i_i_C(2), 0;];
Ja_transl  = t3;
