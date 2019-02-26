% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP9
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
% Datum: 2019-02-26 22:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRPP9_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_jacobig_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:06
% EndTime: 2019-02-26 22:30:06
% DurationCPUTime: 0.03s
% Computational Cost: add. (9->9), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->15)
t131 = sin(pkin(6));
t135 = sin(qJ(1));
t144 = t135 * t131;
t134 = sin(qJ(2));
t143 = t135 * t134;
t137 = cos(qJ(2));
t142 = t135 * t137;
t138 = cos(qJ(1));
t141 = t138 * t131;
t140 = t138 * t134;
t139 = t138 * t137;
t136 = cos(qJ(3));
t133 = sin(qJ(3));
t132 = cos(pkin(6));
t1 = [0, t144, t132 * t142 + t140 (-t132 * t143 + t139) * t133 - t136 * t144, 0, 0; 0, -t141, -t132 * t139 + t143 (t132 * t140 + t142) * t133 + t136 * t141, 0, 0; 1, t132, -t131 * t137, t131 * t134 * t133 - t132 * t136, 0, 0;];
Jg_rot  = t1;
