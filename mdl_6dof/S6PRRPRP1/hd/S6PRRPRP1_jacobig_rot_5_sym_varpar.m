% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRP1_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_jacobig_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:20
% EndTime: 2019-02-26 20:01:20
% DurationCPUTime: 0.02s
% Computational Cost: add. (15->10), mult. (24->19), div. (0->0), fcn. (40->8), ass. (0->15)
t95 = sin(pkin(10));
t96 = sin(pkin(6));
t105 = t95 * t96;
t97 = cos(pkin(10));
t104 = t97 * t96;
t98 = cos(pkin(6));
t99 = sin(qJ(2));
t103 = t98 * t99;
t100 = cos(qJ(2));
t102 = t95 * t100;
t101 = t97 * t100;
t94 = qJ(3) + pkin(11);
t93 = cos(t94);
t92 = sin(t94);
t1 = [0, t105, t98 * t102 + t97 * t99, 0 (-t95 * t103 + t101) * t92 - t93 * t105, 0; 0, -t104, -t98 * t101 + t95 * t99, 0 (t97 * t103 + t102) * t92 + t93 * t104, 0; 0, t98, -t96 * t100, 0, t96 * t99 * t92 - t98 * t93, 0;];
Jg_rot  = t1;
