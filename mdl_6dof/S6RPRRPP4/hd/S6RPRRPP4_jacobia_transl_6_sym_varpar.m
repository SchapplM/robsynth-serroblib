% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPP4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:58:03
% EndTime: 2019-02-26 20:58:03
% DurationCPUTime: 0.12s
% Computational Cost: add. (193->34), mult. (167->45), div. (0->0), fcn. (193->9), ass. (0->30)
t13 = pkin(9) + qJ(3);
t11 = cos(t13);
t31 = r_i_i_C(2) + qJ(5) + pkin(8);
t9 = sin(t13);
t24 = t31 * t9;
t19 = cos(qJ(4));
t8 = t19 * pkin(4) + pkin(3);
t35 = t24 + t11 * t8 + cos(pkin(9)) * pkin(2) + pkin(1);
t14 = qJ(4) + pkin(10);
t10 = sin(t14);
t12 = cos(t14);
t25 = r_i_i_C(3) + qJ(6);
t33 = pkin(5) + r_i_i_C(1);
t34 = t25 * t10 + t33 * t12 + t8;
t17 = sin(qJ(4));
t32 = pkin(4) * t17;
t30 = t11 * t17;
t18 = sin(qJ(1));
t29 = t18 * t10;
t28 = t18 * t12;
t20 = cos(qJ(1));
t27 = t20 * t10;
t26 = t20 * t12;
t22 = pkin(7) + qJ(2) + t32;
t21 = t31 * t11 - t34 * t9;
t4 = t11 * t26 + t29;
t3 = t11 * t27 - t28;
t2 = t11 * t28 - t27;
t1 = t11 * t29 + t26;
t5 = [-t25 * t1 - t35 * t18 - t33 * t2 + t22 * t20, t18, t21 * t20, t25 * t4 - t33 * t3 + (t18 * t19 - t20 * t30) * pkin(4), t20 * t9, t3; t22 * t18 + t35 * t20 + t25 * t3 + t33 * t4, -t20, t21 * t18, t25 * t2 - t33 * t1 + (-t18 * t30 - t19 * t20) * pkin(4), t18 * t9, t1; 0, 0, t34 * t11 + t24 (-t33 * t10 + t25 * t12 - t32) * t9, -t11, t9 * t10;];
Ja_transl  = t5;
