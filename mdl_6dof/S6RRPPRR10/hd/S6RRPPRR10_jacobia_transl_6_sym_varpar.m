% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR10_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR10_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:33:42
% EndTime: 2019-02-26 21:33:42
% DurationCPUTime: 0.11s
% Computational Cost: add. (159->34), mult. (153->47), div. (0->0), fcn. (171->10), ass. (0->29)
t19 = sin(qJ(2));
t21 = cos(qJ(2));
t26 = pkin(2) + r_i_i_C(3) + pkin(9) + pkin(8) + qJ(4);
t25 = t26 * t21;
t18 = pkin(10) + qJ(5);
t14 = sin(t18);
t27 = qJ(3) + pkin(5) * t14 + sin(pkin(10)) * pkin(4);
t36 = -t27 * t19 - pkin(1) - t25;
t15 = cos(t18);
t31 = pkin(5) * t15;
t34 = pkin(7) + t31 + cos(pkin(10)) * pkin(4) + pkin(3);
t16 = qJ(6) + t18;
t12 = sin(t16);
t13 = cos(t16);
t20 = sin(qJ(1));
t22 = cos(qJ(1));
t28 = t22 * t19;
t5 = -t20 * t12 + t13 * t28;
t6 = t12 * t28 + t20 * t13;
t33 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t29 = t20 * t19;
t7 = t22 * t12 + t13 * t29;
t8 = -t12 * t29 + t22 * t13;
t32 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
t30 = r_i_i_C(1) * t13;
t24 = r_i_i_C(1) * t12 + r_i_i_C(2) * t13 + t27;
t23 = -t26 * t19 + t24 * t21;
t11 = t21 * t12 * r_i_i_C(2);
t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) + t36 * t20 + t34 * t22, t23 * t22, t28, t22 * t21 (-t14 * t20 + t15 * t28) * pkin(5) + t33, t33; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t34 * t20 - t36 * t22, t23 * t20, t29, t20 * t21 (t14 * t22 + t15 * t29) * pkin(5) + t32, t32; 0, t24 * t19 + t25, -t21, t19, t11 + (-t30 - t31) * t21, -t21 * t30 + t11;];
Ja_transl  = t1;
