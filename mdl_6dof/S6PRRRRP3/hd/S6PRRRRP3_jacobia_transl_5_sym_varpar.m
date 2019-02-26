% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRP3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:16:20
% EndTime: 2019-02-26 20:16:21
% DurationCPUTime: 0.17s
% Computational Cost: add. (197->41), mult. (384->74), div. (0->0), fcn. (482->12), ass. (0->34)
t25 = sin(qJ(3));
t28 = cos(qJ(3));
t21 = qJ(4) + qJ(5);
t19 = sin(t21);
t20 = cos(t21);
t27 = cos(qJ(4));
t33 = t27 * pkin(4) + r_i_i_C(1) * t20 - r_i_i_C(2) * t19 + pkin(3);
t42 = r_i_i_C(3) + pkin(10) + pkin(9);
t46 = t42 * t25 + t33 * t28 + pkin(2);
t22 = sin(pkin(11));
t26 = sin(qJ(2));
t29 = cos(qJ(2));
t37 = cos(pkin(11));
t38 = cos(pkin(6));
t34 = t38 * t37;
t11 = t22 * t26 - t29 * t34;
t12 = t22 * t29 + t26 * t34;
t23 = sin(pkin(6));
t35 = t23 * t37;
t8 = t12 * t28 - t25 * t35;
t45 = (t11 * t20 - t8 * t19) * r_i_i_C(1) + (-t11 * t19 - t8 * t20) * r_i_i_C(2);
t36 = t22 * t38;
t14 = -t26 * t36 + t37 * t29;
t41 = t23 * t25;
t10 = t14 * t28 + t22 * t41;
t13 = t37 * t26 + t29 * t36;
t44 = (-t10 * t19 + t13 * t20) * r_i_i_C(1) + (-t10 * t20 - t13 * t19) * r_i_i_C(2);
t40 = t23 * t28;
t16 = t38 * t25 + t26 * t40;
t39 = t23 * t29;
t43 = (-t16 * t19 - t20 * t39) * r_i_i_C(1) + (-t16 * t20 + t19 * t39) * r_i_i_C(2);
t24 = sin(qJ(4));
t32 = t24 * pkin(4) + t19 * r_i_i_C(1) + t20 * r_i_i_C(2) + pkin(8);
t1 = [0, -t13 * t46 + t32 * t14, t42 * t10 + t33 * (-t14 * t25 + t22 * t40) (-t10 * t24 + t13 * t27) * pkin(4) + t44, t44, 0; 0, -t11 * t46 + t32 * t12, t42 * t8 + t33 * (-t12 * t25 - t28 * t35) (t11 * t27 - t24 * t8) * pkin(4) + t45, t45, 0; 1 (t32 * t26 + t46 * t29) * t23, t42 * t16 + t33 * (-t26 * t41 + t38 * t28) (-t16 * t24 - t27 * t39) * pkin(4) + t43, t43, 0;];
Ja_transl  = t1;
