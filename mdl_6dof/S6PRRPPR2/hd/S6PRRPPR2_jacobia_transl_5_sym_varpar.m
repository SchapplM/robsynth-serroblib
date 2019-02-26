% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPPR2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:45
% EndTime: 2019-02-26 19:58:46
% DurationCPUTime: 0.12s
% Computational Cost: add. (125->30), mult. (200->52), div. (0->0), fcn. (249->10), ass. (0->29)
t15 = qJ(3) + pkin(11);
t13 = sin(t15);
t14 = cos(t15);
t23 = cos(qJ(3));
t26 = r_i_i_C(3) + qJ(5);
t34 = pkin(4) - r_i_i_C(2);
t35 = pkin(3) * t23 + t26 * t13 + t34 * t14 + pkin(2);
t33 = r_i_i_C(1) + qJ(4) + pkin(8);
t16 = sin(pkin(10));
t17 = sin(pkin(6));
t32 = t16 * t17;
t18 = cos(pkin(10));
t31 = t17 * t18;
t22 = sin(qJ(2));
t30 = t17 * t22;
t29 = t17 * t23;
t19 = cos(pkin(6));
t28 = t19 * t22;
t24 = cos(qJ(2));
t27 = t19 * t24;
t21 = sin(qJ(3));
t10 = -t16 * t28 + t18 * t24;
t9 = t16 * t27 + t18 * t22;
t8 = t16 * t24 + t18 * t28;
t7 = t16 * t22 - t18 * t27;
t5 = t13 * t30 - t14 * t19;
t3 = t10 * t13 - t14 * t32;
t1 = t13 * t8 + t14 * t31;
t2 = [0, t33 * t10 - t35 * t9, t26 * (t10 * t14 + t13 * t32) - t34 * t3 + (-t10 * t21 + t16 * t29) * pkin(3), t9, t3, 0; 0, t33 * t8 - t35 * t7, t26 * (-t13 * t31 + t14 * t8) - t34 * t1 + (-t18 * t29 - t21 * t8) * pkin(3), t7, t1, 0; 1 (t33 * t22 + t35 * t24) * t17, t26 * (t13 * t19 + t14 * t30) - t34 * t5 + (t19 * t23 - t21 * t30) * pkin(3), -t17 * t24, t5, 0;];
Ja_transl  = t2;
