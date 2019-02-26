% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR3
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRR3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:19:52
% EndTime: 2019-02-26 20:19:52
% DurationCPUTime: 0.19s
% Computational Cost: add. (341->50), mult. (495->82), div. (0->0), fcn. (620->14), ass. (0->37)
t32 = sin(qJ(3));
t34 = cos(qJ(3));
t29 = qJ(4) + qJ(5);
t25 = cos(t29);
t19 = pkin(5) * t25 + cos(qJ(4)) * pkin(4);
t26 = qJ(6) + t29;
t22 = sin(t26);
t23 = cos(t26);
t38 = r_i_i_C(1) * t23 - r_i_i_C(2) * t22 + pkin(3) + t19;
t47 = r_i_i_C(3) + pkin(11) + pkin(10) + pkin(9);
t51 = t47 * t32 + t38 * t34 + pkin(2);
t30 = sin(pkin(12));
t33 = sin(qJ(2));
t35 = cos(qJ(2));
t42 = cos(pkin(12));
t43 = cos(pkin(6));
t39 = t43 * t42;
t11 = t30 * t33 - t35 * t39;
t12 = t30 * t35 + t33 * t39;
t31 = sin(pkin(6));
t40 = t31 * t42;
t8 = t12 * t34 - t32 * t40;
t50 = (t11 * t23 - t8 * t22) * r_i_i_C(1) + (-t11 * t22 - t8 * t23) * r_i_i_C(2);
t41 = t30 * t43;
t14 = -t33 * t41 + t42 * t35;
t46 = t31 * t32;
t10 = t14 * t34 + t30 * t46;
t13 = t42 * t33 + t35 * t41;
t49 = (-t10 * t22 + t13 * t23) * r_i_i_C(1) + (-t10 * t23 - t13 * t22) * r_i_i_C(2);
t45 = t31 * t34;
t16 = t43 * t32 + t33 * t45;
t44 = t31 * t35;
t48 = (-t16 * t22 - t23 * t44) * r_i_i_C(1) + (-t16 * t23 + t22 * t44) * r_i_i_C(2);
t24 = sin(t29);
t18 = pkin(5) * t24 + sin(qJ(4)) * pkin(4);
t37 = t22 * r_i_i_C(1) + t23 * r_i_i_C(2) + pkin(8) + t18;
t1 = [0, -t13 * t51 + t37 * t14, t47 * t10 + t38 * (-t14 * t32 + t30 * t45) -t10 * t18 + t13 * t19 + t49 (-t10 * t24 + t13 * t25) * pkin(5) + t49, t49; 0, -t11 * t51 + t37 * t12, t47 * t8 + t38 * (-t12 * t32 - t34 * t40) t11 * t19 - t8 * t18 + t50 (t11 * t25 - t24 * t8) * pkin(5) + t50, t50; 1 (t37 * t33 + t51 * t35) * t31, t47 * t16 + t38 * (-t33 * t46 + t43 * t34) -t16 * t18 - t19 * t44 + t48 (-t16 * t24 - t25 * t44) * pkin(5) + t48, t48;];
Ja_transl  = t1;
