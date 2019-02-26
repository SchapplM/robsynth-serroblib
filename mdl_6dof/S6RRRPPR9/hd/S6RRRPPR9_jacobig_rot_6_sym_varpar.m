% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPPR9_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:10
% EndTime: 2019-02-26 22:08:10
% DurationCPUTime: 0.03s
% Computational Cost: add. (9->9), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->15)
t130 = sin(pkin(6));
t134 = sin(qJ(1));
t143 = t134 * t130;
t133 = sin(qJ(2));
t142 = t134 * t133;
t136 = cos(qJ(2));
t141 = t134 * t136;
t137 = cos(qJ(1));
t140 = t137 * t130;
t139 = t137 * t133;
t138 = t137 * t136;
t135 = cos(qJ(3));
t132 = sin(qJ(3));
t131 = cos(pkin(6));
t1 = [0, t143, t131 * t141 + t139, 0, 0 -(-t131 * t142 + t138) * t132 + t135 * t143; 0, -t140, -t131 * t138 + t142, 0, 0 -(t131 * t139 + t141) * t132 - t135 * t140; 1, t131, -t130 * t136, 0, 0, -t130 * t133 * t132 + t131 * t135;];
Jg_rot  = t1;
