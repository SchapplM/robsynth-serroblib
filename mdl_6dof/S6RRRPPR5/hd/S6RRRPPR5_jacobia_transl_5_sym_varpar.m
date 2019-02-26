% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:05:45
% EndTime: 2019-02-26 22:05:45
% DurationCPUTime: 0.19s
% Computational Cost: add. (226->52), mult. (373->81), div. (0->0), fcn. (469->12), ass. (0->35)
t28 = cos(qJ(3));
t17 = t28 * pkin(3) + pkin(2);
t20 = qJ(3) + pkin(11);
t18 = sin(t20);
t19 = cos(t20);
t21 = sin(pkin(12));
t23 = cos(pkin(12));
t33 = r_i_i_C(1) * t23 - r_i_i_C(2) * t21 + pkin(4);
t39 = r_i_i_C(3) + qJ(5);
t43 = t39 * t18 + t33 * t19 + t17;
t22 = sin(pkin(6));
t26 = sin(qJ(2));
t42 = t22 * t26;
t27 = sin(qJ(1));
t41 = t22 * t27;
t30 = cos(qJ(1));
t40 = t22 * t30;
t38 = cos(pkin(6));
t29 = cos(qJ(2));
t35 = t30 * t38;
t10 = t26 * t35 + t27 * t29;
t37 = t10 * t19 - t18 * t40;
t36 = t27 * t38;
t25 = sin(qJ(3));
t34 = t22 * (pkin(3) * t25 + pkin(8));
t24 = -qJ(4) - pkin(9);
t32 = t21 * r_i_i_C(1) + t23 * r_i_i_C(2) - t24;
t1 = t10 * t18 + t19 * t40;
t12 = -t26 * t36 + t30 * t29;
t11 = t30 * t26 + t29 * t36;
t9 = t27 * t26 - t29 * t35;
t7 = t18 * t42 - t38 * t19;
t6 = t12 * t19 + t18 * t41;
t5 = t12 * t18 - t19 * t41;
t2 = [(-t9 * t21 - t23 * t37) * r_i_i_C(1) + (t21 * t37 - t9 * t23) * r_i_i_C(2) - t37 * pkin(4) - t10 * t17 + t9 * t24 - t27 * pkin(1) - t39 * t1 + t30 * t34, -t11 * t43 + t32 * t12, t39 * t6 + (-t12 * t25 + t28 * t41) * pkin(3) - t33 * t5, t11, t5, 0; (t11 * t21 + t6 * t23) * r_i_i_C(1) + (t11 * t23 - t6 * t21) * r_i_i_C(2) + t6 * pkin(4) + t12 * t17 - t11 * t24 + t30 * pkin(1) + t39 * t5 + t27 * t34, t32 * t10 - t43 * t9, t39 * t37 + (-t10 * t25 - t28 * t40) * pkin(3) - t33 * t1, t9, t1, 0; 0 (t32 * t26 + t43 * t29) * t22, t39 * (t38 * t18 + t19 * t42) + (-t25 * t42 + t38 * t28) * pkin(3) - t33 * t7, -t22 * t29, t7, 0;];
Ja_transl  = t2;
