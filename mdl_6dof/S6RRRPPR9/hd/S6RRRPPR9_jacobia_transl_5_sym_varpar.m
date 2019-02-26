% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR9_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR9_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:09
% EndTime: 2019-02-26 22:08:09
% DurationCPUTime: 0.20s
% Computational Cost: add. (198->54), mult. (504->91), div. (0->0), fcn. (646->10), ass. (0->36)
t31 = sin(qJ(2));
t32 = sin(qJ(1));
t34 = cos(qJ(2));
t35 = cos(qJ(1));
t41 = cos(pkin(6));
t38 = t35 * t41;
t20 = t31 * t38 + t32 * t34;
t30 = sin(qJ(3));
t33 = cos(qJ(3));
t28 = sin(pkin(6));
t46 = t28 * t35;
t10 = t20 * t33 - t30 * t46;
t19 = t32 * t31 - t34 * t38;
t27 = sin(pkin(11));
t29 = cos(pkin(11));
t53 = t10 * t27 - t19 * t29;
t43 = r_i_i_C(2) + qJ(4);
t52 = pkin(3) * t33 + t43 * t30 + pkin(2);
t51 = pkin(4) + r_i_i_C(1);
t49 = t27 * t33;
t48 = t28 * t32;
t47 = t28 * t33;
t45 = t29 * t33;
t44 = t33 * t34;
t42 = r_i_i_C(3) + qJ(5);
t39 = t32 * t41;
t9 = t20 * t30 + t33 * t46;
t36 = -t42 * t27 - t51 * t29 - pkin(3);
t22 = -t31 * t39 + t35 * t34;
t21 = t35 * t31 + t34 * t39;
t18 = t41 * t30 + t31 * t47;
t17 = t28 * t31 * t30 - t41 * t33;
t14 = t22 * t33 + t30 * t48;
t13 = t22 * t30 - t32 * t47;
t3 = t14 * t27 - t21 * t29;
t1 = [pkin(8) * t46 - t32 * pkin(1) - t20 * pkin(2) - t10 * pkin(3) - t19 * pkin(9) + t51 * (-t10 * t29 - t19 * t27) - t43 * t9 - t42 * t53, t22 * pkin(9) + t51 * (-t21 * t45 + t22 * t27) + t42 * (-t21 * t49 - t22 * t29) - t52 * t21, t36 * t13 + t43 * t14, t13, t3, 0; pkin(8) * t48 + t35 * pkin(1) + t22 * pkin(2) + t14 * pkin(3) + t21 * pkin(9) + t51 * (t14 * t29 + t21 * t27) + t42 * t3 + t43 * t13, t20 * pkin(9) + t51 * (-t19 * t45 + t20 * t27) + t42 * (-t19 * t49 - t20 * t29) - t52 * t19, t43 * t10 + t36 * t9, t9, t53, 0; 0 (t51 * (t27 * t31 + t29 * t44) + t42 * (t27 * t44 - t29 * t31) + pkin(9) * t31 + t52 * t34) * t28, t36 * t17 + t43 * t18, t17, t28 * t34 * t29 + t18 * t27, 0;];
Ja_transl  = t1;
