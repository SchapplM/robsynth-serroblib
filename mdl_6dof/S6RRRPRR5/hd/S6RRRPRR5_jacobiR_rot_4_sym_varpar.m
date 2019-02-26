% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR5_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_jacobiR_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:22
% EndTime: 2019-02-26 22:18:22
% DurationCPUTime: 0.02s
% Computational Cost: add. (20->5), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
t76 = cos(qJ(1));
t75 = sin(qJ(1));
t74 = qJ(2) + qJ(3);
t73 = cos(t74);
t72 = sin(t74);
t71 = t76 * t73;
t70 = t76 * t72;
t69 = t75 * t73;
t68 = t75 * t72;
t1 = [t76, 0, 0, 0, 0, 0; t75, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t69, t70, t70, 0, 0, 0; -t71, t68, t68, 0, 0, 0; 0, -t73, -t73, 0, 0, 0; -t68, t71, t71, 0, 0, 0; t70, t69, t69, 0, 0, 0; 0, t72, t72, 0, 0, 0;];
JR_rot  = t1;
