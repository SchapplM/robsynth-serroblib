% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRRRR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobia_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:42
% EndTime: 2019-02-26 19:42:42
% DurationCPUTime: 0.15s
% Computational Cost: add. (198->41), mult. (459->77), div. (0->0), fcn. (599->14), ass. (0->40)
t24 = sin(pkin(13));
t25 = sin(pkin(12));
t28 = cos(pkin(13));
t29 = cos(pkin(12));
t31 = cos(pkin(6));
t40 = t29 * t31;
t16 = -t24 * t25 + t28 * t40;
t26 = sin(pkin(7));
t30 = cos(pkin(7));
t27 = sin(pkin(6));
t41 = t29 * t27;
t13 = -t16 * t26 - t30 * t41;
t23 = qJ(4) + qJ(5);
t21 = sin(t23);
t22 = cos(t23);
t17 = t24 * t40 + t25 * t28;
t33 = sin(qJ(3));
t35 = cos(qJ(3));
t38 = t16 * t30 - t26 * t41;
t8 = t17 * t35 + t38 * t33;
t50 = (t13 * t22 - t21 * t8) * r_i_i_C(1) + (-t13 * t21 - t22 * t8) * r_i_i_C(2);
t45 = t25 * t31;
t19 = -t24 * t45 + t28 * t29;
t18 = -t24 * t29 - t28 * t45;
t44 = t26 * t27;
t37 = t18 * t30 + t25 * t44;
t10 = t19 * t35 + t37 * t33;
t46 = t25 * t27;
t14 = -t18 * t26 + t30 * t46;
t49 = (-t10 * t21 + t14 * t22) * r_i_i_C(1) + (-t10 * t22 - t14 * t21) * r_i_i_C(2);
t42 = t28 * t30;
t43 = t26 * t31;
t12 = t33 * t43 + (t24 * t35 + t33 * t42) * t27;
t15 = -t28 * t44 + t30 * t31;
t48 = (-t12 * t21 + t15 * t22) * r_i_i_C(1) + (-t12 * t22 - t15 * t21) * r_i_i_C(2);
t47 = r_i_i_C(3) + pkin(10) + pkin(9);
t34 = cos(qJ(4));
t39 = pkin(4) * t34 + r_i_i_C(1) * t22 - r_i_i_C(2) * t21 + pkin(3);
t32 = sin(qJ(4));
t1 = [0, t46, t47 * t10 + t39 * (-t19 * t33 + t37 * t35) (-t10 * t32 + t14 * t34) * pkin(4) + t49, t49, 0; 0, -t41, t47 * t8 + t39 * (-t17 * t33 + t38 * t35) (t13 * t34 - t32 * t8) * pkin(4) + t50, t50, 0; 1, t31, t47 * t12 + t39 * (t35 * t43 + (-t24 * t33 + t35 * t42) * t27) (-t12 * t32 + t15 * t34) * pkin(4) + t48, t48, 0;];
Ja_transl  = t1;
