% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPPRR3_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobiR_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:50
% EndTime: 2019-02-26 19:45:50
% DurationCPUTime: 0.03s
% Computational Cost: add. (14->10), mult. (40->24), div. (0->0), fcn. (60->8), ass. (0->15)
t44 = cos(pkin(6));
t45 = sin(qJ(2));
t48 = t44 * t45;
t46 = cos(qJ(2));
t47 = t44 * t46;
t43 = cos(pkin(10));
t42 = cos(pkin(11));
t41 = sin(pkin(6));
t40 = sin(pkin(10));
t39 = sin(pkin(11));
t38 = -t40 * t48 + t43 * t46;
t37 = -t40 * t47 - t43 * t45;
t36 = t40 * t46 + t43 * t48;
t35 = -t40 * t45 + t43 * t47;
t1 = [0, t37 * t42 + t38 * t39, 0, 0, 0, 0; 0, t35 * t42 + t36 * t39, 0, 0, 0, 0; 0 (t39 * t45 + t42 * t46) * t41, 0, 0, 0, 0; 0, -t37 * t39 + t38 * t42, 0, 0, 0, 0; 0, -t35 * t39 + t36 * t42, 0, 0, 0, 0; 0 (-t39 * t46 + t42 * t45) * t41, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
