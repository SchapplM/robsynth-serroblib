% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR11_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR11_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:19
% EndTime: 2019-02-26 21:34:19
% DurationCPUTime: 0.17s
% Computational Cost: add. (259->50), mult. (442->80), div. (0->0), fcn. (561->12), ass. (0->36)
t46 = r_i_i_C(3) + pkin(10);
t45 = pkin(2) + pkin(9) + qJ(4);
t24 = sin(pkin(6));
t27 = sin(qJ(2));
t44 = t24 * t27;
t28 = sin(qJ(1));
t43 = t24 * t28;
t30 = cos(qJ(2));
t42 = t24 * t30;
t31 = cos(qJ(1));
t41 = t24 * t31;
t40 = cos(pkin(6));
t39 = t24 * (pkin(8) + cos(pkin(11)) * pkin(4) + pkin(3));
t38 = sin(pkin(11)) * pkin(4) + qJ(3);
t37 = t28 * t40;
t36 = t31 * t40;
t26 = sin(qJ(6));
t29 = cos(qJ(6));
t35 = t29 * r_i_i_C(1) - t26 * r_i_i_C(2) + pkin(5);
t13 = t28 * t27 - t30 * t36;
t22 = pkin(11) + qJ(5);
t20 = sin(t22);
t21 = cos(t22);
t7 = -t13 * t20 + t21 * t41;
t34 = t13 * t21 + t20 * t41;
t33 = t26 * r_i_i_C(1) + t29 * r_i_i_C(2) + t45;
t32 = t35 * t20 - t46 * t21 + t38;
t16 = -t27 * t37 + t31 * t30;
t15 = t31 * t27 + t30 * t37;
t14 = t27 * t36 + t28 * t30;
t12 = -t20 * t42 + t40 * t21;
t4 = t15 * t20 + t21 * t43;
t3 = -t15 * t21 + t20 * t43;
t2 = t16 * t26 + t4 * t29;
t1 = t16 * t29 - t4 * t26;
t5 = [-t28 * pkin(1) - t38 * t13 - t33 * t14 + t31 * t39 + t46 * t34 + t35 * t7, -t33 * t15 + t32 * t16, t15, t16, -t35 * t3 + t46 * t4, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t31 * pkin(1) + t4 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t38 * t15 + t45 * t16 + t28 * t39 + t46 * t3, -t33 * t13 + t32 * t14, t13, t14, t35 * t34 - t46 * t7 (t14 * t29 + t7 * t26) * r_i_i_C(1) + (-t14 * t26 + t7 * t29) * r_i_i_C(2); 0 (t32 * t27 + t33 * t30) * t24, -t42, t44, t46 * t12 + t35 * (-t40 * t20 - t21 * t42) (-t12 * t26 + t29 * t44) * r_i_i_C(1) + (-t12 * t29 - t26 * t44) * r_i_i_C(2);];
Ja_transl  = t5;
