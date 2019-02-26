% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP9_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP9_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:50:27
% EndTime: 2019-02-26 21:50:28
% DurationCPUTime: 0.20s
% Computational Cost: add. (297->52), mult. (479->82), div. (0->0), fcn. (604->12), ass. (0->39)
t49 = pkin(5) + r_i_i_C(1);
t27 = -pkin(9) - qJ(3);
t28 = sin(qJ(5));
t31 = cos(qJ(5));
t35 = t31 * r_i_i_C(2) + t49 * t28 - t27;
t19 = cos(pkin(11)) * pkin(3) + pkin(2);
t23 = pkin(11) + qJ(4);
t21 = sin(t23);
t22 = cos(t23);
t20 = t31 * pkin(5) + pkin(4);
t37 = t31 * r_i_i_C(1) - t28 * r_i_i_C(2) + t20;
t48 = r_i_i_C(3) + qJ(6) + pkin(10);
t50 = t48 * t21 + t37 * t22 + t19;
t25 = sin(pkin(6));
t29 = sin(qJ(2));
t47 = t25 * t29;
t30 = sin(qJ(1));
t46 = t25 * t30;
t32 = cos(qJ(2));
t45 = t25 * t32;
t33 = cos(qJ(1));
t44 = t25 * t33;
t43 = cos(pkin(6));
t40 = t33 * t43;
t12 = t29 * t40 + t30 * t32;
t4 = t12 * t22 - t21 * t44;
t41 = t30 * t43;
t39 = t25 * (pkin(3) * sin(pkin(11)) + pkin(8));
t13 = t33 * t29 + t32 * t41;
t14 = -t29 * t41 + t33 * t32;
t8 = t14 * t22 + t21 * t46;
t1 = t13 * t31 - t8 * t28;
t3 = t12 * t21 + t22 * t44;
t11 = t30 * t29 - t32 * t40;
t10 = t43 * t21 + t22 * t47;
t9 = t21 * t47 - t43 * t22;
t7 = t14 * t21 - t22 * t46;
t2 = t13 * t28 + t8 * t31;
t5 = [-t30 * pkin(1) - t35 * t11 - t12 * t19 - t48 * t3 + t33 * t39 - t37 * t4, -t13 * t50 + t35 * t14, t13, -t37 * t7 + t48 * t8, -t2 * r_i_i_C(2) + t49 * t1, t7; t33 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t14 * t19 + t8 * t20 + t48 * t7 + t30 * t39 + (t28 * pkin(5) - t27) * t13, -t11 * t50 + t35 * t12, t11, -t37 * t3 + t48 * t4 (-t11 * t28 - t4 * t31) * r_i_i_C(2) + t49 * (t11 * t31 - t4 * t28) t3; 0 (t35 * t29 + t50 * t32) * t25, -t45, t48 * t10 - t37 * t9 (-t10 * t31 + t28 * t45) * r_i_i_C(2) + t49 * (-t10 * t28 - t31 * t45) t9;];
Ja_transl  = t5;
