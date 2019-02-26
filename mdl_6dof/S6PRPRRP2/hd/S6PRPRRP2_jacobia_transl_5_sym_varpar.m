% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRP2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRP2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:50:54
% EndTime: 2019-02-26 19:50:54
% DurationCPUTime: 0.18s
% Computational Cost: add. (198->40), mult. (517->75), div. (0->0), fcn. (679->12), ass. (0->33)
t22 = sin(pkin(11));
t29 = sin(qJ(2));
t32 = cos(qJ(2));
t40 = cos(pkin(11));
t35 = -t29 * t22 + t32 * t40;
t28 = sin(qJ(4));
t31 = cos(qJ(4));
t27 = sin(qJ(5));
t30 = cos(qJ(5));
t37 = t30 * r_i_i_C(1) - t27 * r_i_i_C(2) + pkin(4);
t45 = pkin(9) + r_i_i_C(3);
t33 = t28 * t45 + t37 * t31 + pkin(3);
t23 = sin(pkin(10));
t24 = sin(pkin(6));
t44 = t23 * t24;
t25 = cos(pkin(10));
t43 = t25 * t24;
t26 = cos(pkin(6));
t42 = t26 * t32;
t19 = -t32 * t22 - t29 * t40;
t17 = t19 * t26;
t7 = -t25 * t17 + t23 * t35;
t38 = -t23 * t17 - t25 * t35;
t36 = t27 * r_i_i_C(1) + t30 * r_i_i_C(2) + pkin(8);
t34 = t35 * t26;
t16 = t19 * t24;
t15 = t35 * t24;
t12 = -t16 * t31 + t26 * t28;
t9 = t25 * t19 - t23 * t34;
t6 = t23 * t19 + t25 * t34;
t4 = t28 * t44 - t31 * t38;
t2 = -t28 * t43 + t7 * t31;
t1 = [0 (-t23 * t42 - t25 * t29) * pkin(2) - t36 * t38 + t33 * t9, t44, t45 * t4 + t37 * (t28 * t38 + t31 * t44) (-t4 * t27 - t9 * t30) * r_i_i_C(1) + (t9 * t27 - t4 * t30) * r_i_i_C(2), 0; 0 (-t23 * t29 + t25 * t42) * pkin(2) + t36 * t7 + t33 * t6, -t43, t45 * t2 + t37 * (-t7 * t28 - t31 * t43) (-t2 * t27 - t6 * t30) * r_i_i_C(1) + (-t2 * t30 + t6 * t27) * r_i_i_C(2), 0; 1, t24 * t32 * pkin(2) + t33 * t15 - t36 * t16, t26, t45 * t12 + t37 * (t16 * t28 + t26 * t31) (-t12 * t27 - t15 * t30) * r_i_i_C(1) + (-t12 * t30 + t15 * t27) * r_i_i_C(2), 0;];
Ja_transl  = t1;
