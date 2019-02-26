% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPP9_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP9_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:15
% EndTime: 2019-02-26 22:30:16
% DurationCPUTime: 0.19s
% Computational Cost: add. (312->58), mult. (788->99), div. (0->0), fcn. (1016->10), ass. (0->40)
t35 = sin(qJ(2));
t36 = sin(qJ(1));
t39 = cos(qJ(2));
t49 = cos(pkin(6));
t59 = cos(qJ(1));
t43 = t49 * t59;
t26 = t35 * t43 + t36 * t39;
t34 = sin(qJ(3));
t38 = cos(qJ(3));
t32 = sin(pkin(6));
t45 = t32 * t59;
t14 = t26 * t38 - t34 * t45;
t25 = t36 * t35 - t39 * t43;
t33 = sin(qJ(4));
t37 = cos(qJ(4));
t1 = t14 * t33 - t25 * t37;
t2 = t14 * t37 + t25 * t33;
t48 = pkin(5) + pkin(10) + r_i_i_C(1);
t60 = pkin(3) * t38 + t48 * t34 + pkin(2);
t47 = pkin(4) + r_i_i_C(3) + qJ(6);
t50 = r_i_i_C(2) + qJ(5);
t40 = t50 * t33 + t47 * t37 + pkin(3);
t56 = t32 * t36;
t55 = t32 * t38;
t54 = t32 * t39;
t53 = t33 * t38;
t52 = t37 * t38;
t51 = t38 * t39;
t44 = t36 * t49;
t42 = -t26 * t34 - t38 * t45;
t28 = -t35 * t44 + t59 * t39;
t27 = t59 * t35 + t39 * t44;
t24 = t49 * t34 + t35 * t55;
t18 = t28 * t38 + t34 * t56;
t17 = t28 * t34 - t36 * t55;
t12 = t24 * t37 - t33 * t54;
t11 = t24 * t33 + t37 * t54;
t6 = t18 * t37 + t27 * t33;
t5 = t18 * t33 - t27 * t37;
t3 = [-t36 * pkin(1) - t26 * pkin(2) - t14 * pkin(3) + pkin(8) * t45 - t25 * pkin(9) - t50 * t1 - t47 * t2 + t48 * t42, t28 * pkin(9) + t50 * (-t27 * t53 - t28 * t37) + t47 * (-t27 * t52 + t28 * t33) - t60 * t27, -t40 * t17 + t48 * t18, -t47 * t5 + t50 * t6, t5, t6; t59 * pkin(1) + t28 * pkin(2) + t18 * pkin(3) + pkin(8) * t56 + t27 * pkin(9) + t48 * t17 + t47 * t6 + t50 * t5, t26 * pkin(9) + t50 * (-t25 * t53 - t26 * t37) + t47 * (-t25 * t52 + t26 * t33) - t60 * t25, t48 * t14 + t40 * t42, -t47 * t1 + t50 * t2, t1, t2; 0 (t50 * (t33 * t51 - t35 * t37) + t47 * (t33 * t35 + t37 * t51) + pkin(9) * t35 + t60 * t39) * t32, t48 * t24 + t40 * (-t32 * t35 * t34 + t49 * t38) -t47 * t11 + t50 * t12, t11, t12;];
Ja_transl  = t3;
