% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRP1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:39:38
% EndTime: 2019-02-26 22:39:38
% DurationCPUTime: 0.14s
% Computational Cost: add. (206->33), mult. (162->43), div. (0->0), fcn. (170->10), ass. (0->33)
t24 = qJ(2) + qJ(3);
t21 = qJ(4) + t24;
t18 = sin(t21);
t19 = cos(t21);
t25 = sin(qJ(5));
t44 = r_i_i_C(2) * t25;
t51 = pkin(10) + r_i_i_C(3);
t52 = t18 * t44 + t19 * t51;
t49 = t19 * pkin(4) + t51 * t18;
t27 = cos(qJ(5));
t45 = r_i_i_C(1) * t27;
t33 = (-pkin(4) - t45) * t18;
t30 = t33 - pkin(3) * sin(t24);
t17 = pkin(3) * cos(t24);
t22 = cos(qJ(2)) * pkin(2);
t48 = pkin(1) + t17 + t22 + t49;
t28 = cos(qJ(1));
t41 = t25 * t28;
t26 = sin(qJ(1));
t40 = t26 * t25;
t39 = t26 * t27;
t38 = t27 * t28;
t37 = t52 * t26;
t35 = t52 * t28;
t32 = -sin(qJ(2)) * pkin(2) + t30;
t31 = (-t44 + t45) * t19 + t49;
t29 = t17 + t31;
t23 = -pkin(9) - pkin(8) - pkin(7);
t4 = t19 * t38 + t40;
t3 = -t19 * t41 + t39;
t2 = -t19 * t39 + t41;
t1 = t19 * t40 + t38;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t23 * t28 - t48 * t26, t32 * t28 + t35, t30 * t28 + t35, t28 * t33 + t35, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t26 * t23 + t48 * t28, t32 * t26 + t37, t30 * t26 + t37, t26 * t33 + t37, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, t22 + t29, t29, t31 (-r_i_i_C(1) * t25 - r_i_i_C(2) * t27) * t18, 0;];
Ja_transl  = t5;
