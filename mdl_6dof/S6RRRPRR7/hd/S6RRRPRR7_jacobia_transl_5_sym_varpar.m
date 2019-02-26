% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR7_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR7_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:17
% EndTime: 2019-02-26 22:19:18
% DurationCPUTime: 0.12s
% Computational Cost: add. (204->42), mult. (238->61), div. (0->0), fcn. (289->12), ass. (0->34)
t26 = cos(pkin(6));
t27 = sin(qJ(2));
t30 = cos(qJ(1));
t33 = t30 * t27;
t28 = sin(qJ(1));
t29 = cos(qJ(2));
t34 = t28 * t29;
t10 = t26 * t33 + t34;
t24 = qJ(3) + pkin(12);
t21 = qJ(5) + t24;
t19 = sin(t21);
t25 = sin(pkin(6));
t36 = t25 * t30;
t13 = t19 * t36;
t20 = cos(t21);
t43 = (-t10 * t19 - t20 * t36) * r_i_i_C(1) + (-t10 * t20 + t13) * r_i_i_C(2);
t32 = t30 * t29;
t35 = t28 * t27;
t12 = -t26 * t35 + t32;
t37 = t25 * t28;
t5 = -t12 * t19 + t20 * t37;
t6 = t12 * t20 + t19 * t37;
t42 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t38 = t25 * t27;
t41 = (-t19 * t38 + t26 * t20) * r_i_i_C(1) + (-t26 * t19 - t20 * t38) * r_i_i_C(2);
t15 = pkin(4) * sin(t24) + sin(qJ(3)) * pkin(3);
t40 = pkin(8) + t15;
t39 = r_i_i_C(3) + pkin(10) + qJ(4) + pkin(9);
t16 = pkin(4) * cos(t24) + cos(qJ(3)) * pkin(3);
t14 = pkin(2) + t16;
t31 = t20 * r_i_i_C(1) - t19 * r_i_i_C(2) + t14;
t11 = t26 * t34 + t33;
t9 = -t26 * t32 + t35;
t1 = [-t28 * pkin(1) + t13 * r_i_i_C(1) - t39 * t9 - t31 * t10 + (r_i_i_C(2) * t20 + t40) * t36, -t31 * t11 + t39 * t12, -t12 * t15 + t16 * t37 + t42, t11, t42, 0; t30 * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t39 * t11 + t12 * t14 + t40 * t37, t39 * t10 - t31 * t9, -t10 * t15 - t16 * t36 + t43, t9, t43, 0; 0 (t39 * t27 + t31 * t29) * t25, -t15 * t38 + t26 * t16 + t41, -t25 * t29, t41, 0;];
Ja_transl  = t1;
