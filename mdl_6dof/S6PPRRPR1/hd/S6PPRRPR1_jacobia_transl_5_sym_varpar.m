% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRPR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRRPR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRPR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobia_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:28
% EndTime: 2019-02-26 19:40:29
% DurationCPUTime: 0.12s
% Computational Cost: add. (235->37), mult. (662->71), div. (0->0), fcn. (872->14), ass. (0->41)
t22 = sin(pkin(11));
t24 = sin(pkin(6));
t49 = t22 * t24;
t29 = cos(pkin(6));
t48 = t22 * t29;
t23 = sin(pkin(7));
t47 = t23 * t24;
t46 = t23 * t29;
t26 = cos(pkin(12));
t28 = cos(pkin(7));
t45 = t26 * t28;
t27 = cos(pkin(11));
t44 = t27 * t24;
t43 = t27 * t29;
t42 = r_i_i_C(3) + qJ(5);
t20 = sin(pkin(13));
t25 = cos(pkin(13));
t41 = r_i_i_C(1) * t25 - r_i_i_C(2) * t20 + pkin(4);
t40 = t20 * r_i_i_C(1) + t25 * r_i_i_C(2) + pkin(9);
t21 = sin(pkin(12));
t16 = -t22 * t21 + t26 * t43;
t39 = -t16 * t23 - t28 * t44;
t38 = t16 * t28 - t23 * t44;
t18 = -t27 * t21 - t26 * t48;
t37 = -t18 * t23 + t28 * t49;
t36 = t18 * t28 + t22 * t47;
t35 = -t26 * t47 + t29 * t28;
t30 = sin(qJ(4));
t32 = cos(qJ(4));
t34 = t42 * t30 + t41 * t32 + pkin(3);
t33 = cos(qJ(3));
t31 = sin(qJ(3));
t19 = -t21 * t48 + t27 * t26;
t17 = t21 * t43 + t22 * t26;
t14 = t31 * t46 + (t21 * t33 + t31 * t45) * t24;
t9 = t14 * t30 - t35 * t32;
t8 = t19 * t33 + t36 * t31;
t6 = t17 * t33 + t38 * t31;
t3 = t8 * t30 - t37 * t32;
t1 = t6 * t30 - t39 * t32;
t2 = [0, t49, t40 * t8 + t34 * (-t19 * t31 + t36 * t33) t42 * (t37 * t30 + t8 * t32) - t41 * t3, t3, 0; 0, -t44, t40 * t6 + t34 * (-t17 * t31 + t38 * t33) t42 * (t39 * t30 + t6 * t32) - t41 * t1, t1, 0; 1, t29, t40 * t14 + t34 * (t33 * t46 + (-t21 * t31 + t33 * t45) * t24) t42 * (t14 * t32 + t35 * t30) - t41 * t9, t9, 0;];
Ja_transl  = t2;
