% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% Ja_transl [3x7]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S7RRRRRRR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobia_transl_5_sym_varpar: qJ has to be [7x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S7RRRRRRR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobia_transl_5_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:54:20
% EndTime: 2019-02-26 22:54:20
% DurationCPUTime: 0.24s
% Computational Cost: add. (141->54), mult. (389->104), div. (0->0), fcn. (489->10), ass. (0->47)
t30 = cos(qJ(4));
t24 = sin(qJ(5));
t29 = cos(qJ(5));
t38 = t29 * r_i_i_C(1) - t24 * r_i_i_C(2);
t25 = sin(qJ(4));
t56 = pkin(3) + r_i_i_C(3);
t39 = t56 * t25;
t60 = t38 * t30 + t39;
t26 = sin(qJ(3));
t31 = cos(qJ(3));
t33 = cos(qJ(1));
t42 = t33 * t31;
t28 = sin(qJ(1));
t32 = cos(qJ(2));
t46 = t28 * t32;
t15 = t26 * t46 - t42;
t43 = t33 * t26;
t16 = t31 * t46 + t43;
t27 = sin(qJ(2));
t50 = t27 * t28;
t4 = t16 * t30 + t25 * t50;
t59 = t15 * t29 + t4 * t24;
t58 = t15 * t24 - t4 * t29;
t53 = t26 * t27;
t52 = t26 * t30;
t51 = t26 * t32;
t49 = t27 * t30;
t48 = t27 * t31;
t47 = t27 * t33;
t45 = t32 * t25;
t44 = t32 * t30;
t41 = t24 * t53;
t40 = t29 * t53;
t37 = -t24 * r_i_i_C(1) - t29 * r_i_i_C(2);
t36 = -t16 * t25 + t28 * t49;
t14 = t30 * t48 - t45;
t35 = t25 * t48 + t44;
t34 = t56 * t35;
t20 = -t28 * t26 + t32 * t42;
t19 = -t28 * t31 - t32 * t43;
t18 = t27 * t25 + t31 * t44;
t10 = t14 * t28;
t8 = t20 * t30 + t25 * t47;
t7 = t20 * t25 - t30 * t47;
t2 = t19 * t24 + t8 * t29;
t1 = t19 * t29 - t8 * t24;
t3 = [pkin(2) * t50 + t58 * r_i_i_C(1) + t59 * r_i_i_C(2) + t56 * t36 (-t32 * pkin(2) + t41 * r_i_i_C(1) + t40 * r_i_i_C(2) - t38 * t14 - t34) * t33, t60 * t19 + t37 * t20, -t38 * t7 + t56 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; -pkin(2) * t47 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t56 * t7 (-t10 * t29 + t28 * t41) * r_i_i_C(1) + (t10 * t24 + t28 * t40) * r_i_i_C(2) - pkin(2) * t46 - t28 * t34, -t60 * t15 + t37 * t16, t38 * t36 + t56 * t4, -t59 * r_i_i_C(1) + t58 * r_i_i_C(2), 0, 0; 0 (t18 * t29 - t24 * t51) * r_i_i_C(1) + (-t18 * t24 - t29 * t51) * r_i_i_C(2) - t27 * pkin(2) + t56 * (t31 * t45 - t49) ((-t24 * t31 - t29 * t52) * r_i_i_C(1) + (t24 * t52 - t29 * t31) * r_i_i_C(2) - t26 * t39) * t27, t56 * t14 - t38 * t35 (-t14 * t24 - t40) * r_i_i_C(1) + (-t14 * t29 + t41) * r_i_i_C(2), 0, 0;];
Ja_transl  = t3;
