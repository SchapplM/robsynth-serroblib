% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:31:42
% EndTime: 2019-02-26 22:31:42
% DurationCPUTime: 0.16s
% Computational Cost: add. (248->33), mult. (185->44), div. (0->0), fcn. (196->10), ass. (0->32)
t27 = sin(qJ(6));
t29 = cos(qJ(6));
t55 = r_i_i_C(1) * t27 + r_i_i_C(2) * t29;
t26 = qJ(2) + qJ(3);
t23 = qJ(4) + t26;
t20 = sin(t23);
t21 = cos(t23);
t41 = pkin(4) + pkin(10) + r_i_i_C(3);
t54 = t20 * qJ(5) + t41 * t21;
t52 = (qJ(5) + t55) * t21;
t35 = t41 * t20;
t31 = -t35 - pkin(3) * sin(t26);
t19 = pkin(3) * cos(t26);
t24 = cos(qJ(2)) * pkin(2);
t51 = t19 + t24 + pkin(1) + t54;
t47 = pkin(5) + pkin(9) + pkin(8) + pkin(7);
t28 = sin(qJ(1));
t46 = t28 * t27;
t45 = t28 * t29;
t30 = cos(qJ(1));
t44 = t30 * t27;
t43 = t30 * t29;
t40 = t52 * t28;
t39 = t52 * t30;
t34 = -sin(qJ(2)) * pkin(2) + t31;
t33 = t55 * t20 + t54;
t32 = t19 + t33;
t4 = -t20 * t46 + t43;
t3 = t20 * t45 + t44;
t2 = t20 * t44 + t45;
t1 = t20 * t43 - t46;
t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t51 * t28 + t47 * t30, t34 * t30 + t39, t31 * t30 + t39, -t30 * t35 + t39, t30 * t20, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t47 * t28 + t51 * t30, t34 * t28 + t40, t31 * t28 + t40, -t28 * t35 + t40, t28 * t20, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, t24 + t32, t32, t33, -t21 (-r_i_i_C(1) * t29 + r_i_i_C(2) * t27) * t21;];
Ja_transl  = t5;
