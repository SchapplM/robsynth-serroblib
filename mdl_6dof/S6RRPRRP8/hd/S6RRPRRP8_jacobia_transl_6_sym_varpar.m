% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:49:52
% EndTime: 2019-02-26 21:49:53
% DurationCPUTime: 0.15s
% Computational Cost: add. (266->35), mult. (212->48), div. (0->0), fcn. (244->10), ass. (0->32)
t43 = pkin(5) + r_i_i_C(1);
t35 = r_i_i_C(3) + qJ(6);
t23 = pkin(10) + qJ(4);
t20 = cos(t23);
t11 = pkin(4) * t20 + cos(pkin(10)) * pkin(3) + pkin(2);
t26 = cos(qJ(2));
t24 = sin(qJ(2));
t40 = r_i_i_C(2) + pkin(9) + pkin(8) + qJ(3);
t32 = t40 * t24;
t45 = t26 * t11 + pkin(1) + t32;
t21 = qJ(5) + t23;
t17 = sin(t21);
t18 = cos(t21);
t44 = t35 * t17 + t18 * t43 + t11;
t19 = sin(t23);
t42 = pkin(4) * t19;
t41 = pkin(7) + t42 + sin(pkin(10)) * pkin(3);
t25 = sin(qJ(1));
t38 = t25 * t26;
t27 = cos(qJ(1));
t37 = t27 * t17;
t36 = t27 * t18;
t34 = t35 * t18 * t24;
t33 = t43 * t17;
t7 = t17 * t38 + t36;
t8 = t18 * t38 - t37;
t30 = t35 * t8 - t43 * t7;
t10 = t25 * t17 + t26 * t36;
t9 = -t25 * t18 + t26 * t37;
t29 = t35 * t10 - t43 * t9;
t28 = -t44 * t24 + t40 * t26;
t1 = [-t45 * t25 + t41 * t27 - t35 * t7 - t43 * t8, t28 * t27, t27 * t24 (-t19 * t26 * t27 + t20 * t25) * pkin(4) + t29, t29, t9; t43 * t10 + t41 * t25 + t45 * t27 + t35 * t9, t28 * t25, t25 * t24 (-t19 * t38 - t20 * t27) * pkin(4) + t30, t30, t7; 0, t44 * t26 + t32, -t26 (-t33 - t42) * t24 + t34, -t24 * t33 + t34, t24 * t17;];
Ja_transl  = t1;
