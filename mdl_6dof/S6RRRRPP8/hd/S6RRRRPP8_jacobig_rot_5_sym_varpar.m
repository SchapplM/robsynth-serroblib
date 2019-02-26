% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP8
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRPP8_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_jacobig_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:29:34
% EndTime: 2019-02-26 22:29:34
% DurationCPUTime: 0.03s
% Computational Cost: add. (9->9), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->15)
t134 = sin(pkin(6));
t138 = sin(qJ(1));
t147 = t138 * t134;
t137 = sin(qJ(2));
t146 = t138 * t137;
t140 = cos(qJ(2));
t145 = t138 * t140;
t141 = cos(qJ(1));
t144 = t141 * t134;
t143 = t141 * t137;
t142 = t141 * t140;
t139 = cos(qJ(3));
t136 = sin(qJ(3));
t135 = cos(pkin(6));
t1 = [0, t147, t135 * t145 + t143 (-t135 * t146 + t142) * t136 - t139 * t147, 0, 0; 0, -t144, -t135 * t142 + t146 (t135 * t143 + t145) * t136 + t139 * t144, 0, 0; 1, t135, -t134 * t140, t134 * t137 * t136 - t135 * t139, 0, 0;];
Jg_rot  = t1;
