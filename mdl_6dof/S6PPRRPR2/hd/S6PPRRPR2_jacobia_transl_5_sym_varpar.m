% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRRPR2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRPR2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:56
% EndTime: 2019-02-26 19:40:56
% DurationCPUTime: 0.14s
% Computational Cost: add. (188->35), mult. (524->67), div. (0->0), fcn. (690->12), ass. (0->39)
t44 = pkin(4) - r_i_i_C(2);
t43 = pkin(9) + r_i_i_C(1);
t21 = sin(pkin(11));
t23 = sin(pkin(6));
t42 = t21 * t23;
t27 = cos(pkin(6));
t41 = t21 * t27;
t22 = sin(pkin(7));
t40 = t22 * t23;
t39 = t22 * t27;
t24 = cos(pkin(12));
t26 = cos(pkin(7));
t38 = t24 * t26;
t25 = cos(pkin(11));
t37 = t25 * t23;
t36 = t25 * t27;
t35 = r_i_i_C(3) + qJ(5);
t20 = sin(pkin(12));
t16 = -t21 * t20 + t24 * t36;
t34 = t16 * t26 - t22 * t37;
t18 = -t25 * t20 - t24 * t41;
t33 = t18 * t26 + t21 * t40;
t28 = sin(qJ(4));
t30 = cos(qJ(4));
t32 = t35 * t28 + t44 * t30 + pkin(3);
t31 = cos(qJ(3));
t29 = sin(qJ(3));
t19 = -t20 * t41 + t25 * t24;
t17 = t20 * t36 + t21 * t24;
t15 = -t24 * t40 + t27 * t26;
t14 = -t18 * t22 + t26 * t42;
t13 = -t16 * t22 - t26 * t37;
t12 = t29 * t39 + (t20 * t31 + t29 * t38) * t23;
t9 = t12 * t28 - t15 * t30;
t8 = t19 * t31 + t29 * t33;
t6 = t17 * t31 + t29 * t34;
t3 = -t14 * t30 + t8 * t28;
t1 = -t13 * t30 + t6 * t28;
t2 = [0, t42, t43 * t8 + t32 * (-t19 * t29 + t31 * t33) t35 * (t14 * t28 + t8 * t30) - t44 * t3, t3, 0; 0, -t37, t43 * t6 + t32 * (-t17 * t29 + t31 * t34) t35 * (t13 * t28 + t6 * t30) - t44 * t1, t1, 0; 1, t27, t43 * t12 + t32 * (t31 * t39 + (-t20 * t29 + t31 * t38) * t23) -t44 * t9 + t35 * (t12 * t30 + t15 * t28) t9, 0;];
Ja_transl  = t2;
