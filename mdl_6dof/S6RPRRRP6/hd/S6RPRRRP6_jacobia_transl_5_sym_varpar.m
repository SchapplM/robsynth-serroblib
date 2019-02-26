% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRP6_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP6_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:10:49
% EndTime: 2019-02-26 21:10:50
% DurationCPUTime: 0.11s
% Computational Cost: add. (138->31), mult. (126->43), div. (0->0), fcn. (140->9), ass. (0->31)
t20 = cos(qJ(4));
t10 = pkin(4) * t20 + pkin(3);
t15 = pkin(10) + qJ(3);
t12 = cos(t15);
t11 = sin(t15);
t34 = r_i_i_C(3) + pkin(9) + pkin(8);
t28 = t34 * t11;
t38 = t28 + t12 * t10 + cos(pkin(10)) * pkin(2) + pkin(1);
t16 = qJ(4) + qJ(5);
t14 = cos(t16);
t21 = cos(qJ(1));
t29 = t14 * t21;
t13 = sin(t16);
t19 = sin(qJ(1));
t32 = t13 * t19;
t5 = t12 * t32 + t29;
t30 = t14 * t19;
t31 = t13 * t21;
t6 = -t12 * t30 + t31;
t37 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = -t12 * t31 + t30;
t8 = t12 * t29 + t32;
t36 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t18 = sin(qJ(4));
t35 = pkin(4) * t18;
t33 = t12 * t18;
t27 = pkin(7) + qJ(2) + t35;
t25 = -r_i_i_C(1) * t13 - r_i_i_C(2) * t14;
t24 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + t10;
t23 = -t24 * t11 + t34 * t12;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t38 * t19 + t27 * t21, t19, t23 * t21 (t19 * t20 - t21 * t33) * pkin(4) + t36, t36, 0; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t27 * t19 + t38 * t21, -t21, t23 * t19 (-t19 * t33 - t20 * t21) * pkin(4) + t37, t37, 0; 0, 0, t24 * t12 + t28 (t25 - t35) * t11, t25 * t11, 0;];
Ja_transl  = t1;
