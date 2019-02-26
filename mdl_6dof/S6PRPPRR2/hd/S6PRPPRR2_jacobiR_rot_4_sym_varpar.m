% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPPRR2_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobiR_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:22
% EndTime: 2019-02-26 19:45:22
% DurationCPUTime: 0.03s
% Computational Cost: add. (14->6), mult. (40->16), div. (0->0), fcn. (60->8), ass. (0->13)
t59 = sin(pkin(11));
t62 = cos(pkin(11));
t65 = sin(qJ(2));
t66 = cos(qJ(2));
t67 = t66 * t59 + t65 * t62;
t57 = t65 * t59 - t66 * t62;
t64 = cos(pkin(6));
t63 = cos(pkin(10));
t61 = sin(pkin(6));
t60 = sin(pkin(10));
t56 = t67 * t64;
t55 = t57 * t64;
t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, -t60 * t55 + t63 * t67, 0, 0, 0, 0; 0, t63 * t55 + t60 * t67, 0, 0, 0, 0; 0, t57 * t61, 0, 0, 0, 0; 0, -t60 * t56 - t63 * t57, 0, 0, 0, 0; 0, t63 * t56 - t60 * t57, 0, 0, 0, 0; 0, t67 * t61, 0, 0, 0, 0;];
JR_rot  = t1;
