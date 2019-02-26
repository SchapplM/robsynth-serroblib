% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP2
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
% Datum: 2019-02-26 22:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRP2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:40:09
% EndTime: 2019-02-26 22:40:10
% DurationCPUTime: 0.16s
% Computational Cost: add. (306->34), mult. (253->43), div. (0->0), fcn. (274->10), ass. (0->34)
t28 = sin(qJ(5));
t30 = cos(qJ(5));
t39 = r_i_i_C(3) + qJ(6);
t51 = pkin(5) + r_i_i_C(1);
t57 = t39 * t28 + t51 * t30;
t56 = pkin(10) + r_i_i_C(2);
t27 = qJ(2) + qJ(3);
t24 = qJ(4) + t27;
t22 = cos(t24);
t54 = t22 * t56;
t21 = sin(t24);
t53 = t22 * pkin(4) + t56 * t21;
t34 = (-pkin(4) - t57) * t21;
t32 = t34 - pkin(3) * sin(t27);
t18 = pkin(3) * cos(t27);
t25 = cos(qJ(2)) * pkin(2);
t52 = pkin(1) + t18 + t25 + t53;
t29 = sin(qJ(1));
t49 = t29 * t54;
t43 = t29 * t28;
t42 = t29 * t30;
t31 = cos(qJ(1));
t41 = t31 * t28;
t40 = t31 * t30;
t38 = t31 * t54;
t36 = t57 * t22 + t53;
t35 = t18 + t36;
t33 = -sin(qJ(2)) * pkin(2) + t32;
t26 = -pkin(9) - pkin(8) - pkin(7);
t4 = t22 * t40 + t43;
t3 = t22 * t41 - t42;
t2 = t22 * t42 - t41;
t1 = t22 * t43 + t40;
t5 = [-t39 * t1 - t51 * t2 - t31 * t26 - t52 * t29, t33 * t31 + t38, t32 * t31 + t38, t31 * t34 + t38, -t51 * t3 + t39 * t4, t3; -t29 * t26 + t39 * t3 + t52 * t31 + t51 * t4, t33 * t29 + t49, t32 * t29 + t49, t29 * t34 + t49, -t51 * t1 + t39 * t2, t1; 0, t25 + t35, t35, t36 (-t51 * t28 + t39 * t30) * t21, t21 * t28;];
Ja_transl  = t5;
