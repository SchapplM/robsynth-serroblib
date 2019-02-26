% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR8_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR8_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobia_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:17
% EndTime: 2019-02-26 20:08:17
% DurationCPUTime: 0.11s
% Computational Cost: add. (111->35), mult. (307->66), div. (0->0), fcn. (394->10), ass. (0->33)
t42 = pkin(3) - r_i_i_C(2);
t19 = sin(pkin(7));
t20 = sin(pkin(6));
t41 = t19 * t20;
t23 = cos(pkin(6));
t40 = t19 * t23;
t22 = cos(pkin(7));
t24 = sin(qJ(3));
t39 = t22 * t24;
t26 = cos(qJ(3));
t38 = t22 * t26;
t25 = sin(qJ(2));
t37 = t23 * t25;
t27 = cos(qJ(2));
t36 = t23 * t27;
t35 = t24 * t25;
t34 = t24 * t27;
t33 = t25 * t26;
t32 = t26 * t27;
t31 = r_i_i_C(3) + qJ(4);
t30 = (pkin(9) + r_i_i_C(1)) * t19;
t18 = sin(pkin(12));
t21 = cos(pkin(12));
t13 = -t18 * t25 + t21 * t36;
t29 = -t13 * t22 + t21 * t41;
t15 = -t18 * t36 - t21 * t25;
t28 = t15 * t22 + t18 * t41;
t16 = -t18 * t37 + t21 * t27;
t14 = t18 * t27 + t21 * t37;
t9 = -t26 * t40 + (-t22 * t32 + t35) * t20;
t3 = t16 * t24 - t28 * t26;
t1 = t14 * t24 + t29 * t26;
t2 = [0, t15 * pkin(2) + t42 * (t15 * t26 - t16 * t39) + t31 * (t15 * t24 + t16 * t38) + t16 * t30, t31 * (t16 * t26 + t28 * t24) - t42 * t3, t3, 0, 0; 0, t13 * pkin(2) + t42 * (t13 * t26 - t14 * t39) + t31 * (t13 * t24 + t14 * t38) + t14 * t30, t31 * (t14 * t26 - t29 * t24) - t42 * t1, t1, 0, 0; 1 (t42 * (-t22 * t35 + t32) + t31 * (t22 * t33 + t34) + pkin(2) * t27 + t25 * t30) * t20, -t42 * t9 + t31 * (t24 * t40 + (t22 * t34 + t33) * t20) t9, 0, 0;];
Ja_transl  = t2;
