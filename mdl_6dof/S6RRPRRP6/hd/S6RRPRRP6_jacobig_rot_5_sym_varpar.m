% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRP6_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_jacobig_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:50
% EndTime: 2019-02-26 21:48:50
% DurationCPUTime: 0.10s
% Computational Cost: add. (18->10), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->17)
t131 = sin(pkin(6));
t136 = sin(qJ(1));
t143 = t136 * t131;
t139 = cos(qJ(1));
t142 = t139 * t131;
t130 = sin(pkin(11));
t132 = cos(pkin(11));
t135 = sin(qJ(2));
t138 = cos(qJ(2));
t141 = t138 * t130 + t135 * t132;
t140 = t135 * t130 - t138 * t132;
t137 = cos(qJ(4));
t134 = sin(qJ(4));
t133 = cos(pkin(6));
t127 = t141 * t133;
t126 = t140 * t133;
t1 = [0, t143, 0, -t136 * t126 + t139 * t141 (-t136 * t127 - t139 * t140) * t134 - t137 * t143, 0; 0, -t142, 0, t139 * t126 + t136 * t141 (t139 * t127 - t136 * t140) * t134 + t137 * t142, 0; 1, t133, 0, t140 * t131, t141 * t134 * t131 - t133 * t137, 0;];
Jg_rot  = t1;
