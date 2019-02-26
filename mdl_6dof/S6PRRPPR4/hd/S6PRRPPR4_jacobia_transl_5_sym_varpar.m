% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPPR4_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR4_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:59:59
% EndTime: 2019-02-26 20:00:00
% DurationCPUTime: 0.14s
% Computational Cost: add. (130->38), mult. (344->70), div. (0->0), fcn. (438->10), ass. (0->34)
t24 = sin(qJ(3));
t26 = cos(qJ(3));
t37 = r_i_i_C(2) + qJ(4);
t44 = pkin(3) * t26 + t37 * t24 + pkin(2);
t43 = pkin(4) + r_i_i_C(1);
t20 = sin(pkin(11));
t42 = t20 * t26;
t22 = sin(pkin(6));
t41 = t22 * t24;
t40 = t22 * t26;
t23 = cos(pkin(11));
t39 = t23 * t26;
t27 = cos(qJ(2));
t38 = t26 * t27;
t36 = r_i_i_C(3) + qJ(5);
t35 = cos(pkin(6));
t34 = cos(pkin(10));
t21 = sin(pkin(10));
t32 = t21 * t35;
t31 = t22 * t34;
t30 = t35 * t34;
t28 = -t36 * t20 - t43 * t23 - pkin(3);
t25 = sin(qJ(2));
t16 = t35 * t24 + t25 * t40;
t15 = t25 * t41 - t35 * t26;
t14 = -t25 * t32 + t34 * t27;
t13 = t34 * t25 + t27 * t32;
t12 = t21 * t27 + t25 * t30;
t11 = t21 * t25 - t27 * t30;
t8 = t14 * t26 + t21 * t41;
t7 = t14 * t24 - t21 * t40;
t6 = t12 * t26 - t24 * t31;
t5 = t12 * t24 + t26 * t31;
t1 = [0, t14 * pkin(8) + t43 * (-t13 * t39 + t14 * t20) + t36 * (-t13 * t42 - t14 * t23) - t44 * t13, t28 * t7 + t37 * t8, t7, -t13 * t23 + t8 * t20, 0; 0, t12 * pkin(8) + t43 * (-t11 * t39 + t12 * t20) + t36 * (-t11 * t42 - t12 * t23) - t44 * t11, t28 * t5 + t37 * t6, t5, -t11 * t23 + t6 * t20, 0; 1 (t36 * (t20 * t38 - t23 * t25) + t43 * (t20 * t25 + t23 * t38) + pkin(8) * t25 + t44 * t27) * t22, t28 * t15 + t37 * t16, t15, t22 * t27 * t23 + t16 * t20, 0;];
Ja_transl  = t1;
