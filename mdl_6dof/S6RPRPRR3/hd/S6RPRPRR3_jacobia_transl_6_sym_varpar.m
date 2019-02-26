% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:50:14
% EndTime: 2019-02-26 20:50:14
% DurationCPUTime: 0.12s
% Computational Cost: add. (206->35), mult. (135->49), div. (0->0), fcn. (150->12), ass. (0->29)
t22 = cos(qJ(3));
t21 = sin(qJ(3));
t30 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(4);
t26 = t30 * t21;
t19 = pkin(11) + qJ(5);
t15 = cos(t19);
t9 = pkin(5) * t15 + cos(pkin(11)) * pkin(4) + pkin(3);
t35 = t22 * t9 + pkin(2) + t26;
t17 = qJ(6) + t19;
t11 = sin(t17);
t12 = cos(t17);
t20 = qJ(1) + pkin(10);
t16 = cos(t20);
t14 = sin(t20);
t29 = t14 * t22;
t5 = t11 * t29 + t16 * t12;
t6 = t16 * t11 - t12 * t29;
t34 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t28 = t16 * t22;
t7 = -t11 * t28 + t14 * t12;
t8 = t14 * t11 + t12 * t28;
t33 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t13 = sin(t19);
t32 = pkin(5) * t13;
t31 = pkin(7) + t32 + sin(pkin(11)) * pkin(4);
t25 = -r_i_i_C(1) * t11 - r_i_i_C(2) * t12;
t24 = r_i_i_C(1) * t12 - r_i_i_C(2) * t11 + t9;
t23 = -t24 * t21 + t30 * t22;
t1 = [-sin(qJ(1)) * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t31 * t16 - t35 * t14, 0, t23 * t16, t16 * t21 (-t13 * t28 + t14 * t15) * pkin(5) + t33, t33; cos(qJ(1)) * pkin(1) + t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t31 * t14 + t35 * t16, 0, t23 * t14, t14 * t21 (-t13 * t29 - t15 * t16) * pkin(5) + t34, t34; 0, 1, t24 * t22 + t26, -t22 (t25 - t32) * t21, t25 * t21;];
Ja_transl  = t1;
