% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:28
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRP5_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:28:00
% EndTime: 2019-02-22 10:28:00
% DurationCPUTime: 0.04s
% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
t77 = pkin(9) + qJ(3);
t73 = sin(t77);
t78 = sin(qJ(1));
t85 = t78 * t73;
t76 = pkin(10) + qJ(5);
t74 = cos(t76);
t84 = t78 * t74;
t75 = cos(t77);
t83 = t78 * t75;
t79 = cos(qJ(1));
t82 = t79 * t73;
t81 = t79 * t74;
t80 = t79 * t75;
t72 = sin(t76);
t71 = t78 * t72 + t74 * t80;
t70 = -t72 * t80 + t84;
t69 = t79 * t72 - t74 * t83;
t68 = t72 * t83 + t81;
t1 = [t69, 0, -t73 * t81, 0, t70, 0; t71, 0, -t73 * t84, 0, -t68, 0; 0, 0, t75 * t74, 0, -t73 * t72, 0; t68, 0, t72 * t82, 0, -t71, 0; t70, 0, t72 * t85, 0, t69, 0; 0, 0, -t75 * t72, 0, -t73 * t74, 0; -t85, 0, t80, 0, 0, 0; t82, 0, t83, 0, 0, 0; 0, 0, t73, 0, 0, 0;];
JR_rot  = t1;
