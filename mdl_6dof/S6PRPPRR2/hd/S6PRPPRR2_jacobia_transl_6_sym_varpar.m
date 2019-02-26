% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPPRR2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:12
% EndTime: 2019-02-26 19:45:12
% DurationCPUTime: 0.20s
% Computational Cost: add. (215->45), mult. (560->78), div. (0->0), fcn. (739->12), ass. (0->35)
t33 = cos(qJ(2));
t41 = cos(pkin(11));
t24 = sin(pkin(11));
t30 = sin(qJ(2));
t44 = t30 * t24;
t18 = -t33 * t41 + t44;
t25 = sin(pkin(10));
t27 = cos(pkin(10));
t42 = cos(pkin(6));
t38 = t42 * t41;
t40 = t33 * t42;
t43 = -t24 * t40 - t30 * t38;
t11 = -t27 * t18 + t25 * t43;
t51 = t25 * t18 + t27 * t43;
t29 = sin(qJ(5));
t32 = cos(qJ(5));
t28 = sin(qJ(6));
t31 = cos(qJ(6));
t37 = t31 * r_i_i_C(1) - t28 * r_i_i_C(2) + pkin(5);
t49 = pkin(9) + r_i_i_C(3);
t50 = t37 * t29 - t49 * t32 + qJ(4);
t26 = sin(pkin(6));
t47 = t25 * t26;
t45 = t27 * t26;
t36 = t28 * r_i_i_C(1) + t31 * r_i_i_C(2) + pkin(3) + pkin(8);
t19 = -t33 * t24 - t30 * t41;
t35 = t33 * t38 - t42 * t44;
t17 = t19 * t26;
t16 = t18 * t26;
t13 = t16 * t29 + t42 * t32;
t10 = t27 * t19 - t25 * t35;
t7 = t25 * t19 + t27 * t35;
t4 = t7 * t29 + t32 * t45;
t2 = -t10 * t29 + t32 * t47;
t1 = [0 (-t25 * t40 - t27 * t30) * pkin(2) + t36 * t10 + t50 * t11, t47, -t10, t49 * t2 + t37 * (-t10 * t32 - t29 * t47) (t11 * t31 - t2 * t28) * r_i_i_C(1) + (-t11 * t28 - t2 * t31) * r_i_i_C(2); 0 (-t25 * t30 + t27 * t40) * pkin(2) + t36 * t7 - t50 * t51, -t45, -t7, -t49 * t4 + t37 * (t29 * t45 - t7 * t32) (t4 * t28 - t31 * t51) * r_i_i_C(1) + (t28 * t51 + t4 * t31) * r_i_i_C(2); 1, t26 * t33 * pkin(2) - t36 * t16 - t17 * t50, t42, t16, t49 * t13 + t37 * (t16 * t32 - t42 * t29) (-t13 * t28 - t17 * t31) * r_i_i_C(1) + (-t13 * t31 + t17 * t28) * r_i_i_C(2);];
Ja_transl  = t1;
