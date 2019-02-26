% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPPRR4_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:11
% EndTime: 2019-02-26 21:30:11
% DurationCPUTime: 0.04s
% Computational Cost: add. (18->10), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->17)
t132 = sin(pkin(6));
t137 = sin(qJ(1));
t144 = t137 * t132;
t140 = cos(qJ(1));
t143 = t140 * t132;
t131 = sin(pkin(11));
t133 = cos(pkin(11));
t136 = sin(qJ(2));
t139 = cos(qJ(2));
t142 = t139 * t131 + t136 * t133;
t141 = t136 * t131 - t139 * t133;
t138 = cos(qJ(5));
t135 = sin(qJ(5));
t134 = cos(pkin(6));
t128 = t142 * t134;
t127 = t141 * t134;
t1 = [0, t144, 0, 0, -t137 * t128 - t140 * t141, t135 * t144 - (-t137 * t127 + t140 * t142) * t138; 0, -t143, 0, 0, t140 * t128 - t137 * t141, -t135 * t143 - (t140 * t127 + t137 * t142) * t138; 1, t134, 0, 0, t142 * t132, -t141 * t138 * t132 + t134 * t135;];
Jg_rot  = t1;
