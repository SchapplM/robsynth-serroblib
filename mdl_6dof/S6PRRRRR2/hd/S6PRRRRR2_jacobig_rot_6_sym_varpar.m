% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRRR2_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:19:15
% EndTime: 2019-02-26 20:19:15
% DurationCPUTime: 0.03s
% Computational Cost: add. (29->11), mult. (46->20), div. (0->0), fcn. (75->8), ass. (0->20)
t136 = sin(pkin(12));
t137 = sin(pkin(6));
t146 = t136 * t137;
t141 = cos(qJ(2));
t145 = t137 * t141;
t138 = cos(pkin(12));
t144 = t138 * t137;
t139 = cos(pkin(6));
t140 = sin(qJ(2));
t143 = t139 * t140;
t142 = t139 * t141;
t135 = qJ(3) + qJ(4);
t134 = cos(t135);
t133 = sin(t135);
t132 = t136 * t142 + t138 * t140;
t131 = t136 * t140 - t138 * t142;
t130 = t137 * t140 * t133 - t139 * t134;
t129 = (-t136 * t143 + t138 * t141) * t133 - t134 * t146;
t128 = (t136 * t141 + t138 * t143) * t133 + t134 * t144;
t1 = [0, t146, t132, t132, t129, t129; 0, -t144, t131, t131, t128, t128; 0, t139, -t145, -t145, t130, t130;];
Jg_rot  = t1;
