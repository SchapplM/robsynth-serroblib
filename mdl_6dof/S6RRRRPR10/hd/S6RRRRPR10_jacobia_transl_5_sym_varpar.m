% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR10_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR10_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:48
% EndTime: 2019-02-26 22:35:48
% DurationCPUTime: 0.17s
% Computational Cost: add. (235->40), mult. (359->62), div. (0->0), fcn. (445->10), ass. (0->37)
t52 = pkin(4) - r_i_i_C(2);
t43 = r_i_i_C(3) + qJ(5);
t34 = cos(qJ(3));
t25 = t34 * pkin(3) + pkin(2);
t28 = qJ(3) + qJ(4);
t26 = sin(t28);
t27 = cos(t28);
t53 = t43 * t26 + t52 * t27 + t25;
t51 = r_i_i_C(1) + pkin(10) + pkin(9);
t29 = sin(pkin(6));
t32 = sin(qJ(2));
t50 = t29 * t32;
t33 = sin(qJ(1));
t49 = t29 * t33;
t36 = cos(qJ(1));
t48 = t29 * t36;
t47 = t33 * t32;
t35 = cos(qJ(2));
t46 = t33 * t35;
t45 = t36 * t32;
t44 = t36 * t35;
t30 = cos(pkin(6));
t19 = t30 * t45 + t46;
t8 = t19 * t27 - t26 * t48;
t31 = sin(qJ(3));
t42 = t29 * (pkin(3) * t31 + pkin(8));
t7 = t19 * t26 + t27 * t48;
t41 = t43 * t8 - t52 * t7;
t21 = -t30 * t47 + t44;
t11 = t21 * t26 - t27 * t49;
t12 = t21 * t27 + t26 * t49;
t40 = -t52 * t11 + t43 * t12;
t16 = t26 * t50 - t30 * t27;
t39 = t43 * (t30 * t26 + t27 * t50) - t52 * t16;
t20 = t30 * t46 + t45;
t18 = -t30 * t44 + t47;
t1 = [-t33 * pkin(1) - t51 * t18 - t19 * t25 + t36 * t42 - t43 * t7 - t52 * t8, -t20 * t53 + t51 * t21 (-t21 * t31 + t34 * t49) * pkin(3) + t40, t40, t11, 0; t36 * pkin(1) + t43 * t11 + t52 * t12 + t51 * t20 + t21 * t25 + t33 * t42, -t18 * t53 + t51 * t19 (-t19 * t31 - t34 * t48) * pkin(3) + t41, t41, t7, 0; 0 (t51 * t32 + t53 * t35) * t29 (t30 * t34 - t31 * t50) * pkin(3) + t39, t39, t16, 0;];
Ja_transl  = t1;
