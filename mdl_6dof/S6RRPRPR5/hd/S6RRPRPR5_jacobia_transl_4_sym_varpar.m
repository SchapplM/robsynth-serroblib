% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR5_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR5_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_jacobia_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:07
% EndTime: 2019-02-26 21:40:08
% DurationCPUTime: 0.12s
% Computational Cost: add. (128->40), mult. (313->67), div. (0->0), fcn. (405->10), ass. (0->32)
t38 = r_i_i_C(3) + pkin(9);
t27 = cos(qJ(2));
t37 = t27 * pkin(2);
t22 = cos(pkin(6));
t36 = t22 * t27;
t20 = sin(pkin(6));
t25 = sin(qJ(1));
t35 = t25 * t20;
t28 = cos(qJ(1));
t34 = t28 * t20;
t23 = sin(qJ(4));
t26 = cos(qJ(4));
t19 = sin(pkin(11));
t21 = cos(pkin(11));
t24 = sin(qJ(2));
t31 = t27 * t19 + t24 * t21;
t12 = t31 * t22;
t14 = t24 * t19 - t27 * t21;
t5 = t28 * t12 - t25 * t14;
t33 = t23 * t34 - t5 * t26;
t32 = t25 * t12 + t28 * t14;
t30 = t26 * r_i_i_C(1) - t23 * r_i_i_C(2) + pkin(3);
t29 = t5 * t23 + t26 * t34;
t18 = pkin(1) + t37;
t13 = t22 * t24 * pkin(2) + (-pkin(8) - qJ(3)) * t20;
t11 = t14 * t22;
t10 = t31 * t20;
t7 = t25 * t11 - t28 * t31;
t4 = -t28 * t11 - t25 * t31;
t2 = t23 * t35 - t26 * t32;
t1 = t23 * t32 + t26 * t35;
t3 = [-t5 * pkin(3) + t33 * r_i_i_C(1) + t29 * r_i_i_C(2) - t28 * t13 - t25 * t18 + t38 * t4, -t38 * t32 + (-t24 * t28 - t25 * t36) * pkin(2) + t30 * t7, t35, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; -pkin(3) * t32 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t25 * t13 + t28 * t18 - t38 * t7, t38 * t5 + (-t24 * t25 + t28 * t36) * pkin(2) + t30 * t4, -t34, -t29 * r_i_i_C(1) + t33 * r_i_i_C(2), 0, 0; 0, t38 * t10 + (-t14 * t30 + t37) * t20, t22 (-t10 * t23 + t22 * t26) * r_i_i_C(1) + (-t10 * t26 - t22 * t23) * r_i_i_C(2), 0, 0;];
Ja_transl  = t3;
