% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRPR6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:12
% EndTime: 2019-02-26 19:49:12
% DurationCPUTime: 0.12s
% Computational Cost: add. (162->38), mult. (320->65), div. (0->0), fcn. (406->12), ass. (0->32)
t37 = r_i_i_C(3) + pkin(9) + qJ(5);
t20 = sin(pkin(6));
t24 = sin(qJ(4));
t36 = t20 * t24;
t25 = sin(qJ(2));
t35 = t20 * t25;
t26 = cos(qJ(4));
t34 = t20 * t26;
t27 = cos(qJ(2));
t33 = t20 * t27;
t22 = cos(pkin(6));
t32 = t22 * t25;
t31 = t22 * t27;
t17 = pkin(11) + qJ(6);
t15 = sin(t17);
t16 = cos(t17);
t30 = t16 * r_i_i_C(1) - t15 * r_i_i_C(2) + cos(pkin(11)) * pkin(5) + pkin(4);
t29 = sin(pkin(11)) * pkin(5) + t15 * r_i_i_C(1) + t16 * r_i_i_C(2) + pkin(2) + pkin(8);
t28 = t30 * t24 - t37 * t26 + qJ(3);
t21 = cos(pkin(10));
t19 = sin(pkin(10));
t12 = t22 * t26 - t24 * t33;
t11 = t22 * t24 + t26 * t33;
t10 = -t19 * t32 + t21 * t27;
t9 = t19 * t31 + t21 * t25;
t8 = t19 * t27 + t21 * t32;
t7 = t19 * t25 - t21 * t31;
t4 = t21 * t34 - t7 * t24;
t3 = t21 * t36 + t7 * t26;
t2 = t19 * t34 + t9 * t24;
t1 = t19 * t36 - t9 * t26;
t5 = [0, t28 * t10 - t29 * t9, t9, -t30 * t1 + t37 * t2, t1 (t10 * t16 - t2 * t15) * r_i_i_C(1) + (-t10 * t15 - t2 * t16) * r_i_i_C(2); 0, t28 * t8 - t29 * t7, t7, t30 * t3 - t37 * t4, -t3 (t4 * t15 + t8 * t16) * r_i_i_C(1) + (-t8 * t15 + t4 * t16) * r_i_i_C(2); 1 (t28 * t25 + t29 * t27) * t20, -t33, -t30 * t11 + t37 * t12, t11 (-t12 * t15 + t16 * t35) * r_i_i_C(1) + (-t12 * t16 - t15 * t35) * r_i_i_C(2);];
Ja_transl  = t5;
