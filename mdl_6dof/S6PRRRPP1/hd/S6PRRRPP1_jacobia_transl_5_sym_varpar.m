% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPP1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPP1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:44
% EndTime: 2019-02-26 20:08:44
% DurationCPUTime: 0.19s
% Computational Cost: add. (167->41), mult. (337->74), div. (0->0), fcn. (423->12), ass. (0->34)
t22 = sin(qJ(3));
t25 = cos(qJ(3));
t17 = qJ(4) + pkin(11);
t15 = sin(t17);
t16 = cos(t17);
t24 = cos(qJ(4));
t29 = t24 * pkin(4) + t16 * r_i_i_C(1) - t15 * r_i_i_C(2) + pkin(3);
t38 = r_i_i_C(3) + qJ(5) + pkin(9);
t39 = t38 * t22 + t29 * t25 + pkin(2);
t19 = sin(pkin(6));
t37 = t19 * t22;
t36 = t19 * t25;
t26 = cos(qJ(2));
t35 = t19 * t26;
t34 = cos(pkin(6));
t33 = cos(pkin(10));
t18 = sin(pkin(10));
t32 = t18 * t34;
t31 = t19 * t33;
t30 = t34 * t33;
t21 = sin(qJ(4));
t28 = t21 * pkin(4) + t15 * r_i_i_C(1) + t16 * r_i_i_C(2) + pkin(8);
t23 = sin(qJ(2));
t10 = t34 * t22 + t23 * t36;
t9 = t23 * t37 - t34 * t25;
t8 = -t23 * t32 + t33 * t26;
t7 = t33 * t23 + t26 * t32;
t6 = t18 * t26 + t23 * t30;
t5 = t18 * t23 - t26 * t30;
t4 = t18 * t37 + t8 * t25;
t3 = -t18 * t36 + t8 * t22;
t2 = -t22 * t31 + t6 * t25;
t1 = t6 * t22 + t25 * t31;
t11 = [0, t28 * t8 - t39 * t7, -t29 * t3 + t38 * t4 (-t4 * t15 + t7 * t16) * r_i_i_C(1) + (-t7 * t15 - t4 * t16) * r_i_i_C(2) + (-t4 * t21 + t7 * t24) * pkin(4), t3, 0; 0, t28 * t6 - t39 * t5, -t29 * t1 + t38 * t2 (-t2 * t15 + t5 * t16) * r_i_i_C(1) + (-t5 * t15 - t2 * t16) * r_i_i_C(2) + (-t2 * t21 + t5 * t24) * pkin(4), t1, 0; 1 (t28 * t23 + t39 * t26) * t19, t38 * t10 - t29 * t9 (-t10 * t15 - t16 * t35) * r_i_i_C(1) + (-t10 * t16 + t15 * t35) * r_i_i_C(2) + (-t10 * t21 - t24 * t35) * pkin(4), t9, 0;];
Ja_transl  = t11;
