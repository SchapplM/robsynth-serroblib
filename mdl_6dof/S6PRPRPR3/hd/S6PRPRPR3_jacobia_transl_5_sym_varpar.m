% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRPR3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:43
% EndTime: 2019-02-26 19:47:44
% DurationCPUTime: 0.11s
% Computational Cost: add. (132->29), mult. (339->52), div. (0->0), fcn. (444->10), ass. (0->28)
t38 = pkin(4) - r_i_i_C(2);
t37 = pkin(8) + r_i_i_C(1);
t21 = sin(pkin(10));
t22 = sin(pkin(6));
t36 = t21 * t22;
t24 = cos(pkin(10));
t35 = t24 * t22;
t25 = cos(pkin(6));
t29 = cos(qJ(2));
t34 = t25 * t29;
t33 = r_i_i_C(3) + qJ(5);
t20 = sin(pkin(11));
t23 = cos(pkin(11));
t27 = sin(qJ(2));
t31 = t29 * t20 + t27 * t23;
t16 = t31 * t25;
t17 = t27 * t20 - t29 * t23;
t7 = t24 * t16 - t21 * t17;
t32 = t21 * t16 + t24 * t17;
t26 = sin(qJ(4));
t28 = cos(qJ(4));
t30 = t33 * t26 + t38 * t28 + pkin(3);
t15 = t17 * t25;
t14 = t31 * t22;
t11 = t14 * t26 - t25 * t28;
t3 = -t26 * t32 - t28 * t36;
t1 = t7 * t26 + t28 * t35;
t2 = [0, -t37 * t32 + (-t21 * t34 - t24 * t27) * pkin(2) + t30 * (t21 * t15 - t24 * t31) t36, t33 * (t26 * t36 - t28 * t32) - t38 * t3, t3, 0; 0, t37 * t7 + (-t21 * t27 + t24 * t34) * pkin(2) + t30 * (-t24 * t15 - t21 * t31) -t35, t33 * (-t26 * t35 + t7 * t28) - t38 * t1, t1, 0; 1, t37 * t14 + (pkin(2) * t29 - t17 * t30) * t22, t25, t33 * (t14 * t28 + t25 * t26) - t38 * t11, t11, 0;];
Ja_transl  = t2;
