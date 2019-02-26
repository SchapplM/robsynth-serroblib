% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR9_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:58:46
% EndTime: 2019-02-26 21:58:46
% DurationCPUTime: 0.03s
% Computational Cost: add. (24->11), mult. (31->18), div. (0->0), fcn. (52->8), ass. (0->19)
t132 = sin(pkin(6));
t136 = cos(qJ(2));
t144 = t132 * t136;
t135 = sin(qJ(1));
t143 = t135 * t132;
t134 = sin(qJ(2));
t142 = t135 * t134;
t141 = t135 * t136;
t137 = cos(qJ(1));
t140 = t137 * t132;
t139 = t137 * t134;
t138 = t137 * t136;
t133 = cos(pkin(6));
t131 = pkin(12) + qJ(4) + qJ(5);
t130 = cos(t131);
t129 = sin(t131);
t128 = t133 * t141 + t139;
t127 = -t133 * t138 + t142;
t1 = [0, t143, 0, t128, t128 (-t133 * t142 + t138) * t129 - t130 * t143; 0, -t140, 0, t127, t127 (t133 * t139 + t141) * t129 + t130 * t140; 1, t133, 0, -t144, -t144, t132 * t134 * t129 - t133 * t130;];
Jg_rot  = t1;
