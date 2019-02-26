% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPPR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:08
% EndTime: 2019-02-26 19:58:08
% DurationCPUTime: 0.13s
% Computational Cost: add. (152->33), mult. (253->56), div. (0->0), fcn. (317->12), ass. (0->31)
t17 = qJ(3) + pkin(11);
t15 = sin(t17);
t16 = cos(t17);
t26 = cos(qJ(3));
t18 = sin(pkin(12));
t21 = cos(pkin(12));
t30 = r_i_i_C(1) * t21 - r_i_i_C(2) * t18 + pkin(4);
t34 = r_i_i_C(3) + qJ(5);
t39 = t26 * pkin(3) + t34 * t15 + t30 * t16 + pkin(2);
t19 = sin(pkin(10));
t20 = sin(pkin(6));
t38 = t19 * t20;
t22 = cos(pkin(10));
t37 = t20 * t22;
t25 = sin(qJ(2));
t36 = t20 * t25;
t35 = t20 * t26;
t33 = cos(pkin(6));
t32 = t25 * t33;
t27 = cos(qJ(2));
t31 = t27 * t33;
t29 = t18 * r_i_i_C(1) + t21 * r_i_i_C(2) + pkin(8) + qJ(4);
t24 = sin(qJ(3));
t10 = -t19 * t32 + t22 * t27;
t9 = t19 * t31 + t22 * t25;
t8 = t19 * t27 + t22 * t32;
t7 = t19 * t25 - t22 * t31;
t5 = t15 * t36 - t33 * t16;
t3 = t10 * t15 - t16 * t38;
t1 = t8 * t15 + t16 * t37;
t2 = [0, t29 * t10 - t39 * t9, t34 * (t10 * t16 + t15 * t38) + (-t10 * t24 + t19 * t35) * pkin(3) - t30 * t3, t9, t3, 0; 0, t29 * t8 - t39 * t7, t34 * (-t15 * t37 + t8 * t16) + (-t22 * t35 - t24 * t8) * pkin(3) - t30 * t1, t7, t1, 0; 1 (t29 * t25 + t39 * t27) * t20, t34 * (t33 * t15 + t16 * t36) + (-t24 * t36 + t33 * t26) * pkin(3) - t30 * t5, -t20 * t27, t5, 0;];
Ja_transl  = t2;
