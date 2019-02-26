% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function Jg_rot = S6RRRRPP9_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_jacobig_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:16
% EndTime: 2019-02-26 22:30:16
% DurationCPUTime: 0.03s
% Computational Cost: add. (9->9), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->15)
t128 = sin(pkin(6));
t132 = sin(qJ(1));
t141 = t132 * t128;
t131 = sin(qJ(2));
t140 = t132 * t131;
t134 = cos(qJ(2));
t139 = t132 * t134;
t135 = cos(qJ(1));
t138 = t135 * t128;
t137 = t135 * t131;
t136 = t135 * t134;
t133 = cos(qJ(3));
t130 = sin(qJ(3));
t129 = cos(pkin(6));
t1 = [0, t141, t129 * t139 + t137 (-t129 * t140 + t136) * t130 - t133 * t141, 0, 0; 0, -t138, -t129 * t136 + t140 (t129 * t137 + t139) * t130 + t133 * t138, 0, 0; 1, t129, -t128 * t134, t128 * t131 * t130 - t129 * t133, 0, 0;];
Jg_rot  = t1;
