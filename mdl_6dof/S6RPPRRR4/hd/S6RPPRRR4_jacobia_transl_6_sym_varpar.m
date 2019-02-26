% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRR4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:36:28
% EndTime: 2019-02-26 20:36:28
% DurationCPUTime: 0.13s
% Computational Cost: add. (150->36), mult. (230->50), div. (0->0), fcn. (294->10), ass. (0->32)
t21 = cos(qJ(4));
t20 = cos(qJ(5));
t14 = t20 * pkin(5) + pkin(4);
t17 = qJ(5) + qJ(6);
t15 = sin(t17);
t16 = cos(t17);
t24 = r_i_i_C(1) * t16 - r_i_i_C(2) * t15 + t14;
t19 = sin(qJ(4));
t35 = r_i_i_C(3) + pkin(9) + pkin(8);
t25 = t35 * t19;
t44 = -t24 * t21 - t25;
t43 = t15 * r_i_i_C(1) + t16 * r_i_i_C(2);
t42 = pkin(1) + pkin(2);
t30 = t16 * t21;
t31 = t15 * t21;
t27 = sin(pkin(10));
t28 = cos(pkin(10));
t32 = sin(qJ(1));
t33 = cos(qJ(1));
t7 = -t32 * t27 - t33 * t28;
t8 = t33 * t27 - t32 * t28;
t40 = (-t7 * t16 + t8 * t31) * r_i_i_C(1) + (t7 * t15 + t8 * t30) * r_i_i_C(2);
t5 = t8 * t16 + t7 * t31;
t6 = t8 * t15 - t7 * t30;
t39 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t18 = sin(qJ(5));
t38 = pkin(5) * t18;
t34 = t43 * t19;
t29 = t18 * t21;
t26 = pkin(7) + t38;
t23 = t24 * t19 - t35 * t21;
t1 = [t33 * qJ(2) + (t26 + t43) * t7 + (pkin(3) - t44) * t8 - t42 * t32, t32, 0, t23 * t7 (t20 * t8 + t7 * t29) * pkin(5) + t39, t39; t32 * qJ(2) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t26 * t8 + (-t21 * t14 - pkin(3) - t25) * t7 + t42 * t33, -t33, 0, t23 * t8 (-t20 * t7 + t8 * t29) * pkin(5) + t40, t40; 0, 0, -1, t44, t19 * t38 + t34, t34;];
Ja_transl  = t1;
